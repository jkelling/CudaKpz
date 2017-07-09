#!/usr/bin/awk -f
BEGIN{
	cnt=0
	TMIN=800000
	SKIP=20
}
{
	if($1 > TMIN) {
		cnt +=1;
		if(cnt >= SKIP) {
			print $0; cnt=0
		}
	}
	else print $0
}
