#!/bin/awk -f

BEGIN {}
{
	if ($9==0 && $10==0 && $11==1){print FILENAME; exit}
}
END{}
