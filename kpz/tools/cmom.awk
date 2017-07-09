#!/usr/bin/awk -f 

# OFS="\t"
function isnum(x){return(x==x+0)}
function cmom(n, hn) {
	ret = hn[n];
	bin = 1;
	for(i = 1; i <= n; ++i)
	{
		bin *= (n + 1 - i)/i;
		ret += ((i%2) ? -bin:bin)*(hn[n-i])*(hn[1]**i);
	}
	return ret;
}

BEGIN{prev=0;}
(isnum($1) && prev != $1){
	OFS="\t"
	OFMT="%.12f"
	# print prev, $1, $2, $3, $4, $6, $7, $8, $9, $10, "ref"
	prev=$1;
	hn[0] = 1
	hn[1] = $3
	hn[2] = $4
	hn[3] = $5
	hn[4] = $6
	hn[5] = $7
	hn[6] = $8
	hn[7] = $9
	print $1, $3, cmom(2, hn), cmom(3, hn), cmom(4, hn), cmom(5, hn), cmom(6, hn), cmom(7, hn)
}
