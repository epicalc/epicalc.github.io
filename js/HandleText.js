//HandleText.js
//Given a string of delimited data, produces a summary of columns
//Requires StatFunctions1.js from the OpenEpi js folder


//-------------------------------
function howMany(pattern,hoststr)
  {
  var str2=hoststr.replace(pattern,"");
  return hoststr.length-str2.length;
  }

 //-------------------------------
 function display(theString)
  {
   document.getElementById("individualData").value=theString;
  }
//-------------------------------
function colNamesFromString(dataString)
{
//Returns a string containing column names
//The first line beginning with an underbar character is taken as column names
//If no line begins with an underbar, then the first line with no digit(s)
//is taken to be the column name line
//If neither situation is found, colNamesFromString is returned as an empty string
 var colNames= dataString.match(/^\s*_.*$/m)[0]     //Line begins with underbar
 alert ("colNames has an underbar colNames="+colNames)
 if(colNames !=null)
    {
      //try removing all lines not containing an initial underbar.
     colNames= colNames.replace(/^\s*_+/m,"" )      //remove initial underbar(s)
    }
 alert ("colNames after removing underbar ="+colNames)


if (colNames==null)
{
 //Replace lines beginning with an asterisk.  These are comments.
 colNames=dataString.replace(/^\s*\*.*$/gm,"")
 
 colNames=dataString.replace(/^\*[\d]+.*/gm,"" );  //Remove all lines containing one or more digits
 //alert("dataString="+ dataString);
// colNames=dataString.match(/^.*[,:;\t][\D]+[,:;\t]?.*$/m)
  colNames=dataString.match(/[\D\s\,:;\t]*$/m)

 //find the first line having no number
}
alert("colNames at 45="+colNames)
return colNames;
}
//-------------------------------

function arrayFromString (dataString)
{
  //Given a string of lines containing delimited data, returns an array of columns
  //with each item in a column representing a line in the original string.
var lineArray=new Array();
//var dataString = dataString.replace(/^\s+|\s+$/g, '') ;  //trim spaces and newlines from either end of the string--Google to see source

var leftQuote=/^\s*["']/       //Regular expression patterns for quotes to left or right of a string
var rightQuote=/["']\s*$/
var firstChar=/^\s*./
//var colNameStr='';
var colNameStr=  colNamesFromString(dataString)

dataString=dataString.replace(/^[\D^,;:]*$/gm, "");     //
//This means REPLACE with a NEWLINE, any part of the string that matches the following:
   // Two NewLines without any numeric digit between them--would include most text and empty lines unless they contain one or more commas
   // OR (|) one or more Newlines at the start (^) of the entire string (an initial blank line).

dataString=dataString.replace(/^\n+/, "");  //Remove any newlines at the beginning of the string
dataString=dataString.replace(/\n\n$/g, '\n'); //Replace double newline at end with single
dataString=dataString.replace(/\n$/g, ""); //Remove final newline
//if (colNameStr=="")
//  { colNameStr="None"}
if (dataString !== "")
  {colNameStr+="\n"}
dataString=colNameStr + dataString;       //Column names before rest of string.
lineArray=dataString.split(/\n+/);  //Splits the entire string into an array of lines based on the Newline character.  Will skip blank lines.

var aRow="";

for (i=0; i<lineArray.length;i++)
   {

    aRow=lineArray[i];

    var arrRowB=new Array();
    var arrRowA=new Array();

    if (howMany(/[,;:\t]/g,aRow) >=1)
      {
      //Line is delimited with ,;: or tab.
       arrRowB=aRow.split(/\s*[,;:\t]\s*/);
      }
    else
      {
       //It must be space delimited or have no delimiters (It might be fixed field)
       arrRowA=aRow.split(/\s+/);
       //For each item, see if it was inside a quotation; if so, put it back with the preceding item
      if (/["']/.test(aRow))
      {
        //Line contains quotes
       inQuote=false;

       for (j=0; j< arrRowA.length; j++)
          {
           if (inQuote==false)
             {
               arrRowB.push(arrRowA[j])
             }
           else
             {
               //must be in quote.  Append to last item
               arrRowB[arrRowB.length-1]+= (" " + arrRowA[j])

             }
           if (leftQuote.test(arrRowA[j])==true )
                 {
                   //left quote
                   inQuote=true;
                 }
           if (rightQuote.test(arrRowA[j]))
                  {
                    //right quote
                    inQuote=false;
                  }
           }
        }
        else
        {
           arrRowB=arrRowA;
        }
      // lineArray[i]=arrRowB;

      }
    lineArray[i]=arrRowB;
   }
                  
var numNonBlankLines=lineArray.length;  //Here just for convenience; This is the number of data lines,
                                   //the first dimension of the line array.

//alert(arrayInfo(lineArray));
//The lineArray array is now ready for whatever you want to do with it.
//You can return it from whatever function we are in at the moment  
//return (arrayInfo(lineArray));
return lineArray;
}

//-------------------------------

function arrayInfo(theArray)
{
  //Returns a string containing information about the array of lines containing column data
var nonBlankLines=theArray.length;
var maxCols=0;
var colSummaryStr="";
var dispStr="There are "+nonBlankLines+" non-blank lines.\n";
dispStr+="The first line is line[0]\n";
dispStr+="line[0] is itself an array, with elements line[0,0], line[0,1], etc.\n";
dispStr+="The number of items is line[0].length, in this case--"+theArray[0].length+"\n";
for (var j=0;j<nonBlankLines;j++)
 {
   if (theArray[j].length>maxCols)
   {
     maxCols=theArray[j].length;
   }
 }
for (var c=0;c<maxCols;c++)
  {
    colSummaryStr+= colSummary(theArray,c)
  }
return colSummaryStr;
}
//-------------------------------


var globalSumArr=Array();
function colSummary(colRowArray, colNum)
{
  //Produces summary of a column in the array of lines
  //arranged as columns by their delimiters
  var n=0;
  var sum=0;
  var sumsqrs=0;
  var mean=0;
  var variance=0;
  var SD=0;
  var cell;
  var strcell;
  var freq=new Array();

  var colSumStr;


 for ( r=1;r<colRowArray.length;r++)
  {
     if (colNum<colRowArray[r].length)
        {
          strcell=colRowArray[r][colNum];
          cell=parseFloat(strcell);
          if (strcell=="")
            { strcell="Missing"}

          if (!isNaN(cell) )
             {
              n+=1;
              sum+= cell;
              sumsqrs+= (cell*cell);
             }
          if (!(freq[strcell]))
             {
               freq[strcell]=1;
             }
          else
             {
               freq[strcell]+=1;
             }
        }
 }
  mean= sum/n;
  variance =  (sumsqrs - ((sum*sum)/n) ) / n;
  SD = Math.sqrt(variance);
  sumArr=Array(9);

  colSumStr="Colnum="+(colNum+1)+"\n";
    sumArr[0]=colNum+1;
  colSumStr+="ColName="+colRowArray[0][colNum]+"\n"
    sumArr[1]= colRowArray[0][colNum]
  colSumStr+="n="+n+ "\n";
    sumArr[2]= n;
  colSumStr+="sum="+sum+"\n";
    sumArr[3]= sum;
  colSumStr+="sumsqrs="+sumsqrs+"\n";
    sumArr[4]= sumsqrs;
  colSumStr+="mean="+fmtSigFig(mean,4)+"\n";
    sumArr[5]= fmtSigFig(mean,4);
  colSumStr+="variance="+fmtSigFig(variance,4)+"\n";
    sumArr[6]=fmtSigFig(variance,4);
  colSumStr+="SD="+fmtSigFig(SD,4)+"\n";
    sumArr[7]= fmtSigFig(SD,4);
  if (freq.length >= 1)
     {
       colSumStr+="Freq:\n";
     }
  var sumOfFreq=0;
  var freqNum=new Array();
  var numIdx=0;
  for (var idx in freq)
       {
         colSumStr+=idx + ": "+freq[idx]+"\n";
         freqNum[numIdx]=idx;
         numIdx+=1;
         sumOfFreq+=freq[idx];
       }
   sumArr[8]=freq;
  //alert (colSumStr);
 // summaryRow=[
 // mygrid.addRow((new Date()).valueOf(),[600,'','','',false,'na',false,''],mygrid.getRowIndex(mygrid.getSelectedId()))
 globalSumArr=sumArr;
 //mygrid.addRow((new Date()).valueOf(),sumArr,mygrid.getRowIndex(mygrid.getSelectedId()))

   //Put data into Array and call calculations

 // var module=document.getElementById('panel1')

  var dataFromBin=new Array();
  var metadata=new Array();
  metadata["junk"]="blank";

  var numdata=new Array();
 //  numdata["E1D1"]=11;
 var numerator       //for proportion
 var denominator    //for proportion
   if ((sum==sumsqrs) && (sumsqrs !== 0))
     {
        numerator=sum;     //values were all zeroes and ones
        denominator=n;
     }
   else
     {
        numerator=freq[freqNum["1"]];
        denominator=sumOfFreq;
     }

   numdata["E1D0"]=denominator;              //denominator in proportion
   numdata["E0D0"]=numerator;


   dataFromBin[0]=metadata;
   dataFromBin[1]=numdata;
//parent.parent.dataFromText=dataFromBin;
//parent.loadData=true;
//parent.useOpenEpiEntry(parent.parent.configCmds1, parent.parent.configCmds2,parent.parent.dataFromText);
 // module.OECalculate(dataFromBin);
// var outObjBin=new Output(dataFromBin)



/*
var EWin= parent.EntryWin;
 parent.OECalculate(dataFromBin)
 return colSumStr;
 */
}

//-------------------------------
function handleChange(textID)
{
  //Produces summary data and places it in a
  //field called summaryData in the current
  //web page
  alert("handleChange in "+textID)
  var element= document.getElementById( textID );
  //alert (arrayFromString(element.value));
  var theArray=arrayFromString(element.value);
  showInGrid(arrayInfo(theArray))
  //element.form.summaryData.value=arrayInfo(theArray);
}
//-------------------------------
function handleSelect(textID)
{
//Does
//http://the-stickman.com/web-development/javascript/finding-selection-cursor-position-in-a-textarea-in-internet-explorer/
    alert("handleSelect in "+textID)
   var element = document.getElementById( textID );

if( document.selection ){
	// The current selection
	var range = document.selection.createRange();
	// We'll use this as a 'dummy'
	var stored_range = range.duplicate();
	// Select all text
	stored_range.moveToElementText( element );
	// Now move 'dummy' end point to end point of original range
	stored_range.setEndPoint( 'EndToEnd', range );
	// Now we can calculate start and end points
	element.selectionStart = stored_range.text.length - range.text.length;
	element.selectionEnd = element.selectionStart + range.text.length;
    // alert("250 in document.selection. element.selectionEnd="+ element.selectionEnd);
}
var theArray=new Array()
var strSelection=element.value.substring(element.selectionStart,element.selectionEnd);
if (strSelection != "")
 {
   if (strSelection.indexOf("\n")>0)
   {
      theArray= arrayFromString(strSelection);
   }
   else
   {
      theArray[0]=strSelection;
   }


//alert("theArray.length="+theArray.length)
if (theArray.length==1 )
     {
       //Handle ColSelection since the selection must be on a single line

       var upThroughSelection=element.value.substring(0,element.selectionEnd);
       var upThroughSelectionArray=arrayFromString(upThroughSelection);

       var x=upThroughSelectionArray.length-1;    //Index of line where selection occurred (first line is 0)
       var y=upThroughSelectionArray[x]  //Last line
       var z=y.length-1;     //Index of the item on that line where selection fell.

       var indexOfColumnSelected=z;
       var colStats=colSummary(arrayFromString(element.value),indexOfColumnSelected);
       // alert (colStats);
    //    element.selection.clear() ;
     }
//alert("start="+element.selectionStart + " end=" + element.selectionEnd)
else
     {
if ((element.selectionEnd-element.selectionStart)>0)
       {
         showInGrid(arrayInfo(theArray));
         //element.form.summaryData.value=arrayInfo(theArray);
       }
//alert(theArray);
     }

}
 }