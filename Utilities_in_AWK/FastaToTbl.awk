#!/usr/bin/awk -f
{
	if (substr($1,1,1)==">")
    		if (NR>1)
                	printf "\n%s\t", substr($0,2,length($0)-1)
    		else 
    			printf "%s\t", substr($0,2,length($0)-1)
	else 
        	printf "%s", $0
}END{printf "\n"}
