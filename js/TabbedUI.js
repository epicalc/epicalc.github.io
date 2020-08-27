// JavaScript Document
//For tabbed user interface
var currentPanel;

function writeToWin(content)
{
  var w = window.open();
  var d = w.document;
  d.open();
  d.write(content);
  d.close();
}

function initPanels(panelInicial)
{
  for (pp = 1; pp<6; pp++)
  {
    if (!(pp==5 && device=="smartphone"))
    {
      if (pp!=panelInicial)
      {
        if(device != "smartphone")
        {
          //just for now
         // document.getElementById('panel'+pp).style.display = 'none';
          //alert("display set to none for panel "+pp+ " in showpanel");
        }
      }
      else
      {
        showPanel(pp)
      }
    }
  }
}
function closeMe()
{

  window.location.href = "../OpenEpiMenu.htm";

}

function showPanelAsWindow(panelNum)
{
  var panelWin;
  var HTMLforDiv = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">';
  HTMLforDiv += '<link rel="stylesheet" href="../css/TabbedUI.css" TYPE="text/css" MEDIA="screen">\n';
  HTMLforDiv += '<link rel="stylesheet" href="../css/TabbedUI-print.css" TYPE="text/css" MEDIA="print">\n';
  //HTMLforDiv+='<meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=2">'
  HTMLforDiv += '<head>'
  HTMLforDiv += '<STYLE type=text/css> '
  HTMLforDiv += 'table{width:290px;background-color:#FFFFFF} tr{width:290px;} td{width=290px;}'
  HTMLforDiv += '</STYLE>'

  HTMLforDiv += '</head>';
  HTMLforDiv += '<body>';
  //The first works in Netscape; the second in IE 6
  //Jan 2007
  var tabWidth = 20;

  HTMLforDiv += '<div id="tab0" class="tab" style="left: 2%;" onClick="closeMe();setTimeout(&quot;self.close()&quot;,1000);"> '
  HTMLforDiv += '<strong>'+t("Menu*")+'</strong>'+'</div>';
  HTMLforDiv += '<div id="tab1" class="tab" style="left: 18%;" onClick = "showPanelAsWindow(1);" onMouseOver="hover(this);" onMouseOut="setState(1)">'
  HTMLforDiv += t('Start*')
  HTMLforDiv += '</div>'
  alert("at 66 in TabbedUI.js")
  //HTMLforDiv+='<div id="tab2" class="tab" style="left: 34%;" onClick = "useOpenEpiEntry();self.close()" onMouseOver="hover(this);" onMouseOut="setState(2);">'
  HTMLforDiv += '<div id="tab2" class="tab" style="left: 34%;" onClick = "showPanel(2);" onMouseOver="hover(this);" onMouseOut="setState(2);">'
  HTMLforDiv += t('Enter')
  HTMLforDiv += '</div>'
  HTMLforDiv += '<div id="tab3" class="tab" style="left: 50%;"showPanel(3);" onMouseOver="hover(this);" onMouseOut="setState(3);">'
  HTMLforDiv += t('Results')
  HTMLforDiv += '</div>'
  HTMLforDiv += '<div id="tab4" class="tab" style="left: 70%;" onClick = "showPanel(4);showExWindow();" onMouseOver="hover(this);" onMouseOut="setState(4)">'
  HTMLforDiv += t('Examples')
  HTMLforDiv += '</div>'
  HTMLforDiv += '<div id="tab5" class="tab" style="left: 70%;" onClick = "showPanel(4);showExWindow();" onMouseOver="hover(this);" onMouseOut="setState(4)">'
  HTMLforDiv += t('Help')
  HTMLforDiv += '</div><br />'
  HTMLforDiv += '<table width="100%" border="0" cellspacing="0" cellpadding="2" bgcolor="#CCCCCC">';
  HTMLforDiv += '  <tr color="#FFFFFF">';
  HTMLforDiv += '<div id="banner" style="background-color:#EEEEEE; color:#FFFFFF;">'

 // if(parent.document.getElementById('panel'+panelNum)!=null)
 // {HTMLforDiv += parent.document.getElementById('panel'+panelNum).innerHTML;}
 // else
  HTMLforDiv += document.getElementById('panel'+panelNum).innerHTML;
  HTMLforDiv += '</body>';
  HTMLforDiv += '</html>';


  if (panelWin && panelWin.open) {panelWin.close();} //smartphone allows only one popup window
  panelWin.document.write(HTMLforDiv);
  panelWin.document.close();
  panelWin.focus();
  currentPanel = panelNum;
}
/*
function showStringInDiv(theString,theDivId)
{
  var pTag = document.createElement("p");
  pTag.innerHTML = theString;
  document.getElementById(theDivId).appendChild(pTag);
}
*/

function addStringToDiv(theString, theDivId)
{
  //  $(theDivId).innerHTML+= theString;    //use prototype to write the HTML code for the table
  jQuery('#'+ theDivId).append(theString); //use jQuery to append the item after the div
}

function insertCSS(cssFile)
{
  cssString = "<style media=\"screen,projection\" type=\"text/css\">" + cssFile+ "</style>";
  //document.getElementsByTagName("head")[0].innerHTML+= cssString;
  jQuery('head').append(cssString);
}

function showPanel(panelNum)
{
  var pp = 0;
  if (currentPanel != null)
  {
    hidePanel();
  }
 // $('panel'+panelNum).show();
  document.getElementById('panel'+panelNum).style.visibility = 'visible';

   jQuery('#panel'+panelNum).show();
   document.getElementById('panel'+panelNum).style.display = 'block';
  if (panelNum==2)
  {
    if ($('panel2').innerHTML.length==0)
    {
      var summary = summaryTable().toString();
      summary='<div id="sumtableDiv">' + summary
      addStringToDiv(summary, "panel2");
      tableFromCommands();
      addStringToDiv( '</div>')  //end of sumtableDiv
      if (currConfig()=="smartphone")
       {
        jQuery('#sumtableDiv').css('left', '2%');
        jQuery('#form1').css('width', '50%');
       }
       else
       {
        jQuery('#sumtableDiv').css('left', '10%');
       }
     }
     //Panel 2 was chosen.  If there is an initData function, call it.  This is
     //used by EpiBin to display sample data, and may be extended later.
     if (typeof(initData) == 'function')
       {
       initData();
       }
    }
  currentPanel = panelNum;
  if (currentPanel==4)
   {
    jQuery('#panel4').css('position', 'absolute');
   }
  setState(panelNum);
}
function hidePanel() {
  //hide visible panel, unhilite tab
  //The following line causes Safari to Crash completely.
  // document.getElementById('panel'+currentPanel).style.display = 'none';
 // removeCalendar();
  //alert("display set to none for panel "+currentPanel+ " in hidepanel")

 jQuery('#panel'+ currentPanel).hide();

  document.getElementById('panel'+currentPanel).style.visibility = 'hidden';
  document.getElementById('tab'+currentPanel).style.backgroundColor = '#ffffff';
  document.getElementById('tab'+currentPanel).style.color = 'navy';
}

function setState(tabNum)
{

  if(tabNum==currentPanel)
  {
    document.getElementById('tab'+tabNum).style.backgroundColor = '#666666'; // was'#ddddff'
    document.getElementById('tab'+tabNum).style.color = 'white';
    document.getElementById('tab'+tabNum).style.borderBottomColor = '#666666';
  }
  else
  {
    document.getElementById('tab'+tabNum).style.backgroundColor = '#ffffff';
    document.getElementById('tab'+tabNum).style.color = 'navy';
  }

}

function hover(tab)
{
  tab.style.backgroundColor = '#ddddff';
  // tab.style.borderBottom='none';
}
//end