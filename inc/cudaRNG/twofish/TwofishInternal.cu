#ifndef CUDA_TWOFISCH_INTERNAL_CU
#define CUDA_TWOFISCH_INTERNAL_CU

union U32 {
	u32 i;
	BYTE b[4];
};

#define RS_MOD 0x14D
#define RHO 0x01010101L

/* 
   gcc is smart enough to convert these to roll instructions.  If you want
   to see for yourself, either do gcc -O3 -S, or change the |'s to +'s and 
   see how slow things get (you lose about 30-50 clocks) :).
*/
#define ROL(x,n) (((x) << ((n) & 0x1F)) | ((x) >> (32-((n) & 0x1F))))
#define ROR(x,n) (((x) >> ((n) & 0x1F)) | ((x) << (32-((n) & 0x1F))))

#if BIG_ENDIAN == 1
#define BSWAP(x) (((ROR(x,8) & 0xFF00FF00) | (ROL(x,8) & 0x00FF00FF)))
#else
#define BSWAP(x) (x)
#endif

#define _b(x, N) (((x) >> (N<<3)) & 0xFF)

/* just casting to byte (instead of masking with 0xFF saves *tons* of clocks 
   (around 50) */
#define b0(x) ((BYTE)(x))
/* this saved 10 clocks */
#define b1(x) ((BYTE)((x) >> 8))
/* use byte cast here saves around 10 clocks */
#define b2(x) (BYTE)((x) >> 16)
/* don't need to mask since all bits are in lower 8 - byte cast here saves
   nothing, but hey, what the hell, it doesn't hurt any */
#define b3(x) (BYTE)((x) >> 24)

/* key-schedule and related */
//RS-Matrix multiply

/* 
   multiply two polynomials represented as u32's, actually called with BYTES,
   but since I'm not really going to too much work to optimize key setup (since
   raw encryption speed is what I'm after), big deal.
*/
inline __device__ u32 polyMult(u32 a, u32 b)
{
    u32 t=0;
    while (a)
    {
		/*printf("A=%X  B=%X  T=%X\n", a, b, t);*/
		if (a&1) t^=b;
		b <<= 1;
		a >>= 1;
    }
    return t;
}
	    
/* take the polynomial t and return the t % modulus in GF(256) */
inline __device__ u32 gfMod(u32 t, u32 modulus)
{
    int i;
    u32 tt;

    modulus <<= 7;
    for (i = 0; i < 8; i++)
    {
		tt = t ^ modulus;
		if (tt < t) t = tt;
		modulus >>= 1;
    }
    return t;
}

/*multiply a and b and return the modulus */
#define gfMult(a, b, modulus) gfMod(polyMult(a, b), modulus)

/* return a u32 containing the result of multiplying the RS Code matrix
   by the sd matrix
*/

inline __device__ u32 RSMatrixMultiply(BYTE sd[8])
{
    int j, k;
    BYTE t;
    U32 result;

    for (j = 0; j < 4; j++)
    {
		t = 0;
		for (k = 0; k < 8; k++)
		{
			/*printf("t=%X  %X\n", t, gfMult(RS[j][k], sd[k], RS_MOD));*/
			t ^= gfMult(RS[j][k], sd[k], RS_MOD);
		}
		result.b[3-j] = t;
    }
    return result.i;
}

#define BYTES_TO_U32(r0, r1, r2, r3) ((r0 << 24) ^ (r1 << 16) ^ (r2 << 8) ^ r3)
#define multEF(x) (x*0xef)
#define mult5B(x) (x*0x5b)

/* the Zero-keyed h function (used by the key setup routine) */
inline __device__ u32 h(u32 X, u32 L[4])
{
    BYTE y0, y1, y2, y3;
    BYTE z0, z1, z2, z3;
    y0 = b0(X);
    y1 = b1(X);
    y2 = b2(X);
    y3 = b3(X);

	    y0 = Q1[y0] ^ b0(L[3]);
	    y1 = Q0[y1] ^ b1(L[3]);
	    y2 = Q0[y2] ^ b2(L[3]);
	    y3 = Q1[y3] ^ b3(L[3]);

	    y0 = Q1[y0] ^ b0(L[2]);
	    y1 = Q1[y1] ^ b1(L[2]);
	    y2 = Q0[y2] ^ b2(L[2]);
	    y3 = Q0[y3] ^ b3(L[2]);

	    y0 = Q1[  Q0 [ Q0[y0] ^ b0(L[1]) ] ^ b0(L[0]) ];
	    y1 = Q0[  Q0 [ Q1[y1] ^ b1(L[1]) ] ^ b1(L[0]) ];
	    y2 = Q1[  Q1 [ Q0[y2] ^ b2(L[1]) ] ^ b2(L[0]) ];
	    y3 = Q0[  Q1 [ Q1[y3] ^ b3(L[1]) ] ^ b3(L[0]) ];

    /* inline the MDS matrix multiply */
    z0 = multEF(y0) ^ y1 ^         multEF(y2) ^ mult5B(y3); 
    z1 = multEF(y0) ^ mult5B(y1) ^ y2 ^         multEF(y3); 
    z2 = mult5B(y0) ^ multEF(y1) ^ multEF(y2) ^ y3; 
    z3 = y0 ^         multEF(y1) ^ mult5B(y2) ^ mult5B(y3); 

    return BYTES_TO_U32(z0, z1, z2, z3);
}

/* given the Sbox keys, create the fully keyed QF */
inline __device__ void fullKey(u32 L[4], u32 QF[4][256])
{
    BYTE y0, y1, y2, y3;

	#ifdef TWOFISH_MULTI
    const int i = threadIdx.x;
	#endif
    
    /* for all input values to the Q permutations */
	#ifndef TWOFISH_MULTI
		if(threadIdx.x == 0)
    		for (int i=0; i<256; ++i)
	#else
		#ifdef TWOFISH_SMALL_BLOCK
			for(int i = threadIdx.x; i < 256; i += blockDim.x)
		#else
    		if(i < 256)
		#endif
	#endif
    {
		/* run the Q permutations */
		y0 = i; y1=i; y2=i; y3=i;

    	y0 = Q1[y0] ^ b0(L[3]);
    	y1 = Q0[y1] ^ b1(L[3]);
    	y2 = Q0[y2] ^ b2(L[3]);
    	y3 = Q1[y3] ^ b3(L[3]);

    	y0 = Q1[y0] ^ b0(L[2]);
    	y1 = Q1[y1] ^ b1(L[2]);
    	y2 = Q0[y2] ^ b2(L[2]);
    	y3 = Q0[y3] ^ b3(L[2]);

    	y0 = Q1[  Q0 [ Q0[y0] ^ b0(L[1]) ] ^ b0(L[0]) ];
    	y1 = Q0[  Q0 [ Q1[y1] ^ b1(L[1]) ] ^ b1(L[0]) ];
    	y2 = Q1[  Q1 [ Q0[y2] ^ b2(L[1]) ] ^ b2(L[0]) ];
    	y3 = Q0[  Q1 [ Q1[y3] ^ b3(L[1]) ] ^ b3(L[0]) ];
	
		/* now do the partial MDS matrix multiplies */
		QF[0][i] = ((multEF(y0) << 24) 
		    | (multEF(y0) << 16) 
		    | (mult5B(y0) << 8)
		    | y0);
		QF[1][i] = ((y1 << 24) 
		    | (mult5B(y1) << 16) 
		    | (multEF(y1) << 8)
		    | multEF(y1));
		QF[2][i] = ((multEF(y2) << 24) 
		    | (y2 << 16) 
		    | (multEF(y2) << 8)
		    | mult5B(y2));
		QF[3][i] = ((mult5B(y3) << 24) 
		    | (multEF(y3) << 16)
		    | (y3 << 8) 
		    | mult5B(y3));
    }
	#ifdef TWOFISH_MULTI
    __syncthreads();
	#endif
}

/* the key schedule routine */
//inline __device__ void keySched(u32 M[], u32 S[4], u32 K[40])
/* key will be copied into Mo[] and Me[] before S[] is written
 * won't be needed afterwards: can be overridden by with S[] */
inline __device__ void keySched(u32 S[], u32 K[40])
{
    __shared__ u32 Mo[4], Me[4];
    __shared__ BYTE vector[8];

	const int k = 4;
	#ifdef TWOFISH_MULTI
	const int i = threadIdx.x;
	#endif
    
	#ifndef TWOFISH_MULTI
    for (int i = 0; i < k; ++i)
	#else
	if(i < k)
	#endif
    {
		Me[i] = BSWAP(S[2*i]);
		Mo[i] = BSWAP(S[2*i+1]);
	}
	#ifdef TWOFISH_MULTI
	__syncthreads();
	#endif
	
	#ifndef TWOFISH_MULTI
	if(threadIdx.x == 0)
    for (int i = 0; i < k; ++i)
	#else
	if(i < 16)
	#endif
    {
		#ifndef TWOFISH_MULTI
		for (int j = 0; j < 4; ++j) vector[j] = _b(Me[i], j);
		for (int j = 0; j < 4; ++j) vector[j+4] = _b(Mo[i], j);
		S[4-i-1] = RSMatrixMultiply(vector);
		#else
		const int j = i&3;
		vector[j] = _b(Me[i>>2], j);
		vector[j+4] = _b(Mo[i>>2], j);
		if(!j)
			S[4-1-(i>>2)] = RSMatrixMultiply(vector);
		#endif
    }
	#ifdef TWOFISH_MULTI
    __syncthreads();
	#endif
	
	#ifndef TWOFISH_MULTI
	if(threadIdx.x == 0)
    for (int i = 0; i < 20; ++i)
	#else
	if(i < 20)
	#endif
    {
		u32 A = h(2*i*RHO, Me);
		u32 B = ROL(h(2*i*RHO + RHO, Mo), 8);
		K[2*i] = A+B;
		K[2*i+1] = ROL(A + 2*B, 9);
    }
}

/* encryption */

/* fully keyed h (aka g) function */
#define fkh(X) (S[0][b0(X)]^S[1][b1(X)]^S[2][b2(X)]^S[3][b3(X)])

/* one encryption round */
#define ENC_ROUND(R0, R1, R2, R3, round, k1,k2) \
    T0 = fkh(R0); \
    T1 = fkh(ROL(R1, 8)); \
    R2 = ROR(R2 ^ (T1 + T0 + k1), 1); \
    R3 = ROL(R3, 1) ^ (2*T1 + T0 + k2); 

__device__ inline void encryptblock(u32 K[40], u32 S[4][256], u32 PT[])
{
    u32 R0, R1, R2, R3;
    u32 T0, T1;

    /* load/byteswap/whiten input */
    R3 = K[3] ^ BSWAP(PT[3]);
    R2 = K[2] ^ BSWAP(PT[2]);
    R1 = K[1] ^ BSWAP(PT[1]);
    R0 = K[0] ^ BSWAP(PT[0]);

    ENC_ROUND(R0, R1, R2, R3, 0, K[8], K[9]);
    ENC_ROUND(R2, R3, R0, R1, 1, K[10], K[11]);
    ENC_ROUND(R0, R1, R2, R3, 2, K[12], K[13]);
    ENC_ROUND(R2, R3, R0, R1, 3, K[14], K[15]);
    ENC_ROUND(R0, R1, R2, R3, 4, K[16], K[17]);
    ENC_ROUND(R2, R3, R0, R1, 5, K[18], K[19]);
    ENC_ROUND(R0, R1, R2, R3, 6, K[20], K[21]);
    ENC_ROUND(R2, R3, R0, R1, 7, K[22], K[23]);
    ENC_ROUND(R0, R1, R2, R3, 8, K[24], K[25]);
    ENC_ROUND(R2, R3, R0, R1, 9, K[26], K[27]);
    ENC_ROUND(R0, R1, R2, R3, 10, K[28], K[29]);
    ENC_ROUND(R2, R3, R0, R1, 11, K[30], K[31]);
    ENC_ROUND(R0, R1, R2, R3, 12, K[32], K[33]);
    ENC_ROUND(R2, R3, R0, R1, 13, K[34], K[35]);
    ENC_ROUND(R0, R1, R2, R3, 14, K[36], K[37]);
    ENC_ROUND(R2, R3, R0, R1, 15, K[38], K[39]);

    /* load/byteswap/whiten output */
    PT[3] = BSWAP(R1 ^ K[7]);
    PT[2] = BSWAP(R0 ^ K[6]);
    PT[1] = BSWAP(R3 ^ K[5]);
    PT[0] = BSWAP(R2 ^ K[4]);
}

/* one decryption round */
#define DEC_ROUND(R0, R1, R2, R3, round) \
    T0 = fkh(R0); \
    T1 = fkh(ROL(R1, 8)); \
    R2 = ROL(R2, 1) ^ (T0 + T1 + K[2*round+8]); \
    R3 = ROR(R3 ^ (T0 + 2*T1 + K[2*round+9]), 1); 

__device__ inline void decryptblock(u32 K[40], u32 S[4][256], u32 PT[])
{
    u32 T0, T1;
    u32 R0, R1, R2, R3;

    /* load/byteswap/whiten input */
    R3 = K[7] ^ BSWAP(PT[3]);
    R2 = K[6] ^ BSWAP(PT[2]);
    R1 = K[5] ^ BSWAP(PT[1]);
    R0 = K[4] ^ BSWAP(PT[0]);

    DEC_ROUND(R0, R1, R2, R3, 15);
    DEC_ROUND(R2, R3, R0, R1, 14);
    DEC_ROUND(R0, R1, R2, R3, 13);
    DEC_ROUND(R2, R3, R0, R1, 12);
    DEC_ROUND(R0, R1, R2, R3, 11);
    DEC_ROUND(R2, R3, R0, R1, 10);
    DEC_ROUND(R0, R1, R2, R3, 9);
    DEC_ROUND(R2, R3, R0, R1, 8);
    DEC_ROUND(R0, R1, R2, R3, 7);
    DEC_ROUND(R2, R3, R0, R1, 6);
    DEC_ROUND(R0, R1, R2, R3, 5);
    DEC_ROUND(R2, R3, R0, R1, 4);
    DEC_ROUND(R0, R1, R2, R3, 3);
    DEC_ROUND(R2, R3, R0, R1, 2);
    DEC_ROUND(R0, R1, R2, R3, 1);
    DEC_ROUND(R2, R3, R0, R1, 0);

    /* load/byteswap/whiten output */
    PT[3] = BSWAP(R1 ^ K[3]);
    PT[2] = BSWAP(R0 ^ K[2]);
    PT[1] = BSWAP(R3 ^ K[1]);
    PT[0] = BSWAP(R2 ^ K[0]);

}

static __device__ void init(u32* dKey, u32 K[40], u32 QF[4][256])
{
	__shared__ u32 S[8];
	
	const int gid = blockIdx.x*8; // 8 = keysize
	if(threadIdx.x < 8)
	{
		S[threadIdx.x] = dKey[gid+threadIdx.x];
	}
	#ifdef TWOFISH_MULTI
	__syncthreads();
	#endif

	#ifndef TWOFISH_MULTI
	if(threadIdx.x == 0)
	#endif
	{
	keySched(S, K);
	#ifdef TWOFISH_MULTI
	__syncthreads();
	#endif
	fullKey(S, QF);
	}
}

#endif
