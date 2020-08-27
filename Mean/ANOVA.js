// ANOVA.JS, a JavaScript module for ANOVA.HTM
//Here is the text that goes on the first page of the application, the one presented by running the HTM.  
var Title="Analysis of Variance (ANOVA)";

var Authors="<b>Statistics</b><br>Minn M. Soe and Kevin M. Sullivan, Emory University<br>"+
"<br>"+
"<b>Interface</b><br>"+
"Andrew G. Dean, EpiInformatics.com, "+
"<br>and Roger A. Mir<br> ";
      
var Description=
"This one-way ANOVA module compares the means of two or more independent samples. Entering sample size, mean and standard deviation (or standard error) of each sample group will test for significant difference among the sample means. The confidence intervals for each individual mean are also displayed.";

//The text in the next variables will be inserted into the HTML document that comes up in response to the Exercise link
var Demo="Suppose a study was conducted to analyse the difference in mean finger-wrist tapping score (MAXFWT) by lead-exposure group. Group 1 is the control group with no history of lead exposure, group 2 the currently exposed group with elevated blood-lead level and group 5 the previously exposed group. To compare the mean scores in above three groups, one-way ANOVA has to be used."+
	  "<ul>"+
	    "<li>To analyse this, enter the given values in respective cells in Open Epi ANOVA for comparing two or more means, and click on Calculate</li>"+
		"<li>In the new window screen, p-value from <i>F</i> test would be 0.012.</li>"+
      "</ul>"+
	  "Conclusion: There is an overall significant difference among mean MAXFWT scores in the three groups."+
	  "<p>Reference: Bernard Rosner. Fundamentals of Biostatistics(5th edition). (Data obtained from table 12.5; pg 533).";
	  

var Exercises="currently not available";
 
//----------------------------------------------------------------------------------------------------;
	  
var Pi=Math.PI; var PiD2=Pi/2; var PiD4=Pi/4; var Pi2=2*Pi;
var e=2.718281828459045235; var e10 = 1.105170918075647625;
var Deg=180/Pi;

//function keys;
function ChiSq(x,n) {
    if(n==1 & x>24) {return "<0.001"};
    
    if(x>1000 | n>100000) {
        var q=ChiSq((x-n)*(x-n)/(2*n),1)/2
        if(x>n) {return q} {return 1-q}
        }
        
    var Pi = 3.141592653589793238462643383;
    var p=Math.exp(-0.5*x); if((n%2)==1) { p=p*Math.sqrt(2*x/Pi) }
    var k=n; while(k>=2) { p=p*x/k; k=k-2 }
    
    var t=p; var a=n; while(t>0.0000000001*p) { a=a+2; t=t*x/a; p=p+t };
    return 1-p;
	}
	
function StudT(t,n) {
    t=Math.abs(t); var w=t/Math.sqrt(n); var th=Math.atan(w)
    if(n==1) { return 1-th/PiD2 }
    var sth=Math.sin(th); var cth=Math.cos(th)
    if((n%2)==1)
        { return 1-(th+sth*cth*StatCom(cth*cth,2,n-3,-1))/PiD2 }
        else
        { return 1-sth*StatCom(cth*cth,1,n-3,-1) }
    }
function StatCom(q,i,j,b) {
    var zz=1; var z=zz; var k=i; while(k<=j) { zz=zz*q*k/(k-b); z=z+zz; k=k+2 }
    return z
    }
function AStudT(p,n) { var v=0.5; var dv=0.5; var t=0
    while(dv>1e-6) { t=1/v-1; dv=dv/2; if(StudT(t,n)>p) { v=v-dv } else { v=v+dv } }
    return t;
    }

function FishF(f,n1,n2) {
    var x=n2/(n1*f+n2)
    if((n1%2)==0) { return StatCom(1-x,n2,n1+n2-4,n2-2)*Math.pow(x,n2/2) }
    if((n2%2)==0){ return 1-StatCom(x,n1,n1+n2-4,n1-2)*Math.pow(1-x,n1/2) }
    var th=Math.atan(Math.sqrt(n1*f/n2)); var a=th/PiD2; var sth=Math.sin(th); var cth=Math.cos(th)
    if(n2>1) { a=a+sth*cth*StatCom(cth*cth,2,n2-3,-1)/PiD2 }
    if(n1==1) { return 1-a }
    var c=4*StatCom(sth*sth,n2+1,n1+n2-4,n2-2)*sth*Math.pow(cth,n2)/Pi
    if(n2==1) { return 1-a+c/2 }
    var k=2; while(k<=(n2-1)/2) {c=c*k/(k-.5); k=k+1 }
    return 1-a+c
    }
function AFishF(p,n1,n2) { var v=0.5; var dv=0.5; var f=0
    while(dv>1e-10) { f=1/v-1; dv=dv/2; if(FishF(f,n1,n2)>p) { v=v-dv } else { v=v+dv } }
    return f
    }

function space(n) {
	if (n==0) {n=""};
	return n;
}

// create arrays;

var n=new Array(11);
var m = new Array(11);
var sd = new Array(11);
var se = new Array(11);
var va = new Array(11);
var l = new Array(11);
var u = new Array(11);
var ll = new Array(11);
var uu = new Array(11);

function CalcAnova(data)
{
//----get data;

var indexname;
for (var i=0;i<10;i++)
  {
  indexname="E"+i+"D0"
  n[i+1] = parseFloat(data[1][indexname]);
  indexname="E"+i+"D1"
  m[i+1] = parseFloat(data[1][indexname]);
  indexname="E"+i+"D2"
  sd[i+1] = parseFloat(data[1][indexname]);
  indexname="E"+i+"D4"
  se[i+1] = parseFloat(data[1][indexname]);
  }



// ------------------------------ANOVA test calculation;
var g=0; var sn=0; var st=0; var sq=0; var sv=0; var sumnum=0; sumnum1=0; sumnum2=0; var sumdenom=0; var cnum=0; 

for (i=1; i<=10; i++) { var ni = n[i]
	if(ni>0){
		if (se[i]!=0) {sd[i] = se[i]*Math.sqrt(n[i]);}  //if SE is entered, it will be converted to SD;
		
		// summing up for Bartlett test;
		va[i]	  = (sd[i]*sd[i]);
		sumnum 	 += (n[i]-1) * (va[i]);
		sumdenom += (n[i]-1);
		
		sumnum1  += (n[i]-1) * Math.log(va[i]);
		cnum 	 += 1/(n[i]-1);
		
		// summing for ANOVA test;
		var ti = ni*m[i];
		var qi = ni*m[i]*m[i] + (ni-1)*sd[i]*sd[i];
		var vi = ti*ti/ni; 	// partA of vSSb;
		
		g = g + 1;			//#of comparison groups;
		sn = sn + ni;
		st = st + ti; 		// partB of vssb;
		sq = sq + qi;		
		sv = sv + vi;
		}
	}

// ANOVA F test;
var vSSb = sv - st*st/sn;  	//SSbet;
var vSSw = sq - sv; 		//SSwithin;
var vSSt = vSSb + vSSw; 	//SStotal;

var vDFb = g - 1;			// degree of freedom;
var vDFt = sn - 1;
var vDFw = vDFt - vDFb;

var vEVb = vSSb / vDFb;		// mean square between;
var vEVw = vSSw / vDFw;		// mean square within; 

// F and p values from ANOVA;
var vF = vEVb / vEVw;
var vP = FishF( vF , vDFb , vDFw );



// --------------------------Bartlett's Homogeneity test;
var vpool= sumnum/sumdenom;
var sumnum2  = (sn-g) * Math.log(vpool);
var cval=  1  +   ( 1/(3*g - 3)) * (cnum - (1/(sn-g)) );

var cs= (sumnum2-sumnum1)/cval; 		//chisquare statistics;
var dfbart=g-1;  				//degree of freedom;
var pbart= ChiSq(cs,dfbart); 			//p value of chisq;

//alert(vF+"  "+vP+"   "+cs+"   "+pbart+" "+dfbart)


// ----------------------------calculating 95% confidence intervals of individual MEANs;
for (i=1; i<=10; i++) {
l[i]=""; u[i]=""; ll[i]=""; uu[i]="";

if(n[i]>0){
if (sd[i]!=0) {se[i] = sd[i]/Math.sqrt(n[i]);}  //if SD is entered, it will be converted to SE to calculate CI;

//based on t-statistics;
var t=AStudT(0.05,n[i]-1); // 0.05 means t statistics at alpha=0.05;
 //Added -1 above Feb 26, 2008
l[i]=m[i] - (t*se[i]); //based on individual Mean;
l[i]=fmtSigFig(l[i],6);

u[i]=m[i] + (t*se[i]);
u[i]=fmtSigFig(u[i],6);

ll[i]=m[i] - t*Math.sqrt(vEVw/n[i]);//based on average mean, assuming equal variance;
ll[i]=fmtSigFig(ll[i],6);

uu[i]=m[i] + t*Math.sqrt(vEVw/n[i]);
uu[i]=fmtSigFig(uu[i],6);
}
}

// Output table;

with (outTable)
{
newtable(8,90);	 //6 columns and 90 pixels per column
title("<h3>" + Title+ "</h3>");
newrow("span8:bold:c:Input Data");
newrow("span8:l: &nbsp;  ");
newrow("color#ffff99:span2:bold:c: Group","color#ffff99:bold:c: N (count)","color#ffff99:bold:c:Mean","color#ffff99:bold:c: Std. Dev.","color#ffff99:bold:c:Std. error");
for (var i=1;i<11;i++)
 {
 if (n[i])
 {newrow("color#ffff99:span2:bold:c:"+i,"c:"+space(n[i]),"c:"+space(m[i]),"c:"+fmtSigFig(sd[i],4),"c:"+fmtSigFig(se[i],4),"bold:l:");}
 }
 newrow();

 newrow("span8:bold:c:ANOVA Table");
 newrow("span8:bold:c: &nbsp;");
 newrow("span2:bold:r:Source of variation","span2:bold:r: Sum of squares","bold:c: d.f","bold:c:Mean square", 	"bold:c:<i>F </i>statistics","bold:c: p-value<tt><sup>1</sup>","color#ffff99:bold:l:");
 newrow("span2:bold:r:Between Groups:","span2:r:"+fmtSigFig(vSSb,6),"c:"+fmtSigFig(vDFb,6),"c:"+fmtSigFig(vEVb,6),"c:"+fmtSigFig(vF,6),"c:"+fmtSigFig(vP,6));
 newrow("span2:bold:r:Within Groups:","span2:r:"+fmtSigFig(vSSw,6),"c:"+fmtSigFig(vDFw,6),"c:"+fmtSigFig(vEVw,6));
 newrow("span2:bold:r:Total:","span2:r:"+fmtSigFig(vSSt,6),"c:"+fmtSigFig(vDFt,6));
 newrow("span8:bold:c: &nbsp;");
	// Bartlett's test result;
 newrow("span3:bold:r:", "bold:r: Chi square","bold:c:d.f","bold:c:p-value<tt><sup>1</sup>","span2:bold:c:");
 newrow("span3:bold:c:Test for equality of variance:","r:"+fmtSigFig(cs,6),"c:"+fmtSigFig(dfbart,6),"c:"+fmtSigFig(pbart,6),"span2:c:");
line(8);
	
	//CI results;
newrow("span2:bold:c:","bold:c:","span2:bold:c:95% CI of individual sample mean","bold:c:","span2:bold:c:95% CI assuming equal variance");
newrow("span8:bold:c: &nbsp;");
newrow("span2:bold:c:Group","bold:l: Mean","bold:c: Lower Limit","bold:c:Upper Limit","bold:c:","bold:c:Lower Limit","bold:c: Upper Limit");
newrow("span2:bold:c:");

for (var i=1;i<11;i++)
 {
 if (m[i])
  {
   newrow("span2:bold:c:"+i,"l:"+m[i],"c:"+l[i],"c:"+u[i],"bold:c:","c:"+ll[i],"c:"+uu[i]);
  }
 }
newrow("span8:bold:c:");
line(8);
newrow("span8:l: <tt><sup>1</sup></tt> <small> p-value (two-tailed)");
endtable();
}
}
//end of calculate ANOVA routine






