INPUT=$1

# INPUT EXAMPLE
# SRR2534126.4587512      112     1       1
# SRR2534126.4587513      124     1       1
# SRR2534126.4587514      91      1       1
# SRR2534126.4587515      92      1       1
# SRR2534126.4587516      86      1       1
# SRR2534126.2655950      102             0
# SRR2534126.4587517      98      1       1
# SRR2534126.2655951      64      4       1
# SRR2534126.4587518      43      1       1
# SRR2534126.2655952      104     1       1


cat $INPUT | awk -v FS="\t" '{

if($2>= 20)
  total++;

# indexes

# 20
# 50
# 100
# 150
# 200
# 250

if($2 >= 20 && $2 < 50){
	if($4!=0){ # it is mapped
		mapped["20"]++;
		if($3==1){
			unique["20"]++;
		}else{
		   not_unique["20"]++;
		}
	}else{
		unmapped["20"]++;
	}
}else if($2 < 100){
	if($4!=0){ 
      mapped["50"]++;
      if($3==1){
			unique["50"]++;
		}else{
		   not_unique["50"]++;
		}
   }else{
     	unmapped["50"]++;
   }
}else if($2 < 150){
	if($4!=0){ 
   	mapped["100"]++;
      if($3==1){
			unique["100"]++;
		}else{
		   not_unique["100"]++;
		}
   }else{
   	unmapped["100"]++;
   }
}else if($2 < 200){
	if($4!=0){ 
   	mapped["150"]++;
      if($3==1){
			unique["150"]++;
		}else{
		   not_unique["150"]++;
		}
        }else{
               	unmapped["150"]++;
        }
}else if($2 < 250){
	if($4!=0){ 
           	mapped["200"]++;
           	if($3==1){
			unique["200"]++;
		}else{
		   not_unique["200"]++;
		}
        }else{
              	unmapped["200"]++;
        }
}else if($2 >= 250){
	if($4!=0){ 
           	mapped["250"]++;
           	if($3==1){
			unique["250"]++;
		}else{
		   not_unique["250"]++;
		}
        }else{
              	unmapped["250"]++;
        }
}
} END {
print "Absolute:"
print "[20,50[",mapped["20"],unique["20"];
print "[50,100[", mapped["50"],unique["50"];
print "[100,150[",mapped["100"],unique["100"];
print "[150,200[",mapped["150"],unique["150"];
print "[200,250[",mapped["200"],unique["200"];
print "[250,...[",mapped["250"],unique["250"];

if(total == mapped["20"] + mapped["50"] + mapped["100"] + mapped["150"] + mapped["200"] + mapped["250"] + unmapped["20"] + unmapped["50"] + unmapped["100"] + unmapped["150"] + unmapped["200"] + unmapped["250"])
	print "Ok.";

print "Relative:"

 
print "[20,50[",mapped["20"]*100/(unmapped["20"]+mapped["20"]), unique["20"]*100/(unmapped["20"]+not_unique["20"] + unique["20"]);
print "[50,100[", mapped["50"]*100/(unmapped["50"]+mapped["50"]), unique["50"]*100/(unmapped["50"]+not_unique["50"] +unique["50"]);
print "[100,150[",mapped["100"]*100/(unmapped["100"]+mapped["100"]), unique["100"]*100/(unmapped["100"]+not_unique["100"] + unique["100"]);
print "[150,200[",mapped["150"]*100/(unmapped["150"]+mapped["150"]), unique["150"]*100/(unmapped["150"]+not_unique["150"] + unique["150"]);
print "[200,250[",mapped["200"]*100/(unmapped["200"]+mapped["200"]), unique["200"]*100/(unmapped["200"]+not_unique["200"] + unique["200"]);
print "[250,...[",mapped["250"]*100/(unmapped["250"]+mapped["250"]), unique["250"]*100/(unmapped["250"]+not_unique["250"] + unique["250"]);


print "Total"

print "[20,50[", (unmapped["20"]+not_unique["20"] + unique["20"])*100/total;
print "[50,100[", (unmapped["50"]+not_unique["50"] +unique["50"])*100/total;
print "[100,150[", (unmapped["100"]+not_unique["100"] + unique["100"])*100/total;
print "[150,200[", (unmapped["150"]+not_unique["150"] + unique["150"])*100/total;
print "[200,250[", (unmapped["200"]+not_unique["200"] + unique["200"])*100/total;
print "[250,...[", (unmapped["250"]+not_unique["250"] + unique["250"])*100/total;

}'


