// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true ];
var arrayMetadata    = [ [ "1", "C4", "0.993694327126503", "C4", "C", "Cold", "Brown" ], [ "2", "C5", "0.638847309833442", "C5", "C", "Heat", "Brown" ], [ "3", "C9", "0.993502086136475", "C9", "C", "Control", "Brown" ], [ "4", "D4", "0.779839904696833", "D4", "D", "Control", "Brown" ], [ "5", "D5", "1.59767277667185", "D5", "D", "Cold", "Brown" ], [ "6", "D6", "1.00225664428161", "D6", "D", "Heat", "Brown" ], [ "7", "E10", "0.993694327126503", "E10", "E", "Cold", "White" ], [ "8", "E11", "0.798507333978348", "E11", "E", "Control", "White" ], [ "9", "E3", "1.26052585366265", "E3", "E", "Heat", "White" ], [ "10", "F7", "0.571928003966705", "F7", "F", "Control", "Brown" ], [ "11", "F8", "0.783488364509574", "F8", "F", "Cold", "Brown" ], [ "12", "F9", "1.57170979901146", "F9", "F", "Heat", "Brown" ], [ "13", "H1", "0.919295418588938", "H1", "X", "Control", "White" ], [ "14", "I1", "2.3106732796896", "I1", "I", "Control", "Brown" ], [ "15", "I2", "1.13469175153054", "I2", "I", "Cold", "Brown" ], [ "16", "I3", "0.724930251597052", "I3", "I", "Heat", "Brown" ], [ "17", "J13", "1.84664123720199", "J13", "J", "Control", "Brown" ], [ "18", "J14", "2.02083237960262", "J14", "J", "Cold", "Brown" ], [ "19", "J15", "0.674692753999186", "J15", "J", "Heat", "Brown" ], [ "20", "K2", "1.11267548244783", "K2", "X", "Cold", "White" ], [ "21", "K3", "0.630028761944678", "K3", "X", "Heat", "White" ], [ "22", "L6", "1.08061072535715", "L6", "Y", "Control", "Brown" ], [ "23", "M2", "0.61867506517791", "M2", "Y", "Cold", "Brown" ], [ "24", "M3", "0.823797190507774", "M3", "Y", "Heat", "Brown" ], [ "25", "N4", "1.33610352763842", "N4", "Z", "Control", "White" ], [ "26", "O2", "1.21530717528284", "O2", "Z", "Cold", "White" ], [ "27", "Q1", "0.801831063521916", "Q1", "Q", "Control", "White" ], [ "28", "Q11", "1.34866714078107", "Q11", "Q", "Heat", "White" ], [ "29", "Q4", "2.97077160101405", "Q4", "Q", "Cold", "White" ], [ "30", "R7", "0.925972571969535", "R7", "R", "Control", "Brown" ], [ "31", "R8", "1.33537750046033", "R8", "R", "Cold", "Brown" ], [ "32", "R9", "0.48207125803266", "R9", "R", "Heat", "Brown" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
