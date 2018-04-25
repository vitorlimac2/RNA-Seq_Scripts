INPUT=$1

cat $INPUT | awk 'NR%4==2{

	l=length($0)

	if(l>=20 && l<50){
		l50++;
	}else if(l<100){
		l100++;
	}else if(l < 150){
		l150++;
	}else if(l < 200){
		l200++;
	}else if(l < 250){
		l250++;
	}else if(l>=250){
		l300++;
	}
}

END {
	total = l50+l100+l150+l200+l250+l300;
	print "[20,50[",l50;
	print "[50,100[",l100;
	print "[100,150[",l150;
	print "[150,200[",l200;
	print "[200,250[",l250;
	print "[250,...]",l300;
	print "Total", total;

	print "###### Density"

	print "[20,50[",l50*100/total;
        print "[50,100[",l100*100/total;
        print "[100,150[",l150*100/total;
        print "[150,200[",l200*100/total;
        print "[200,250[",l250*100/total;
        print "[250,...]",l300*100/total;
        print "Total", total;

}'
