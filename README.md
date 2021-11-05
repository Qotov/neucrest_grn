# neucrest_grn

<html>
<head>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/vis/4.16.1/vis.css" type="text/css" />
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/vis/4.16.1/vis-network.min.js"> </script>
<center>
<h1></h1>
</center>

<!-- <link rel="stylesheet" href="../node_modules/vis/dist/vis.min.css" type="text/css" />
<script type="text/javascript" src="../node_modules/vis/dist/vis.js"> </script>-->

<style type="text/css">

        #mynetwork {
            width: 100%;
            height: 750px;
            background-color: #ffffff;
            border: 1px solid lightgray;
            position: relative;
            float: left;
        }

        
        #loadingBar {
            position:absolute;
            top:0px;
            left:0px;
            width: 100%;
            height: 750px;
            background-color:rgba(200,200,200,0.8);
            -webkit-transition: all 0.5s ease;
            -moz-transition: all 0.5s ease;
            -ms-transition: all 0.5s ease;
            -o-transition: all 0.5s ease;
            transition: all 0.5s ease;
            opacity:1;
        }

        #bar {
            position:absolute;
            top:0px;
            left:0px;
            width:20px;
            height:20px;
            margin:auto auto auto auto;
            border-radius:11px;
            border:2px solid rgba(30,30,30,0.05);
            background: rgb(0, 173, 246); /* Old browsers */
            box-shadow: 2px 0px 4px rgba(0,0,0,0.4);
        }

        #border {
            position:absolute;
            top:10px;
            left:10px;
            width:500px;
            height:23px;
            margin:auto auto auto auto;
            box-shadow: 0px 0px 4px rgba(0,0,0,0.2);
            border-radius:10px;
        }

        #text {
            position:absolute;
            top:8px;
            left:530px;
            width:30px;
            height:50px;
            margin:auto auto auto auto;
            font-size:22px;
            color: #000000;
        }

        div.outerBorder {
            position:relative;
            top:400px;
            width:600px;
            height:44px;
            margin:auto auto auto auto;
            border:8px solid rgba(0,0,0,0.1);
            background: rgb(252,252,252); /* Old browsers */
            background: -moz-linear-gradient(top,  rgba(252,252,252,1) 0%, rgba(237,237,237,1) 100%); /* FF3.6+ */
            background: -webkit-gradient(linear, left top, left bottom, color-stop(0%,rgba(252,252,252,1)), color-stop(100%,rgba(237,237,237,1))); /* Chrome,Safari4+ */
            background: -webkit-linear-gradient(top,  rgba(252,252,252,1) 0%,rgba(237,237,237,1) 100%); /* Chrome10+,Safari5.1+ */
            background: -o-linear-gradient(top,  rgba(252,252,252,1) 0%,rgba(237,237,237,1) 100%); /* Opera 11.10+ */
            background: -ms-linear-gradient(top,  rgba(252,252,252,1) 0%,rgba(237,237,237,1) 100%); /* IE10+ */
            background: linear-gradient(to bottom,  rgba(252,252,252,1) 0%,rgba(237,237,237,1) 100%); /* W3C */
            filter: progid:DXImageTransform.Microsoft.gradient( startColorstr='#fcfcfc', endColorstr='#ededed',GradientType=0 ); /* IE6-9 */
            border-radius:72px;
            box-shadow: 0px 0px 10px rgba(0,0,0,0.2);
        }
        

        

        
        /* position absolute is important and the container has to be relative or absolute as well. */
	    div.popup {
            position:absolute;
            top:0px;
            left:0px;
            display:none;
            background-color:#f5f4ed;
            -moz-border-radius: 3px;
            -webkit-border-radius: 3px;
            border-radius: 3px;
            border: 1px solid #808074;
            box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);
	    }

	    /* hide the original tooltip */
	    .vis-network-tooltip {
	      display:none;
	    }
        
</style>

</head>

<body>
<div id = "mynetwork"></div>

<div id="loadingBar">
    <div class="outerBorder">
        <div id="text">0%</div>
        <div id="border">
            <div id="bar"></div>
        </div>
    </div>
</div>


<script type="text/javascript">

    // initialize global variables.
    var edges;
    var nodes;
    var network; 
    var container;
    var options, data;

    
    // This method is responsible for drawing the graph, returns the drawn network
    function drawGraph() {
        var container = document.getElementById('mynetwork');
        
        

        // parsing and collecting nodes and edges from the python
        nodes = new vis.DataSet([{"id": "TFCP2L1", "label": "TFCP2L1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TFCP2L1\"\u003eTFCP2L1 function\u003c/a\u003e"}, {"id": "SLC16A3", "label": "SLC16A3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC16A3\"\u003eSLC16A3 function\u003c/a\u003e"}, {"id": "GRHL3", "label": "GRHL3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=GRHL3\"\u003eGRHL3 function\u003c/a\u003e"}, {"id": "ATP6AP2", "label": "ATP6AP2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ATP6AP2\"\u003eATP6AP2 function\u003c/a\u003e"}, {"id": "BSG", "label": "BSG", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=BSG\"\u003eBSG function\u003c/a\u003e"}, {"id": "CLDN4", "label": "CLDN4", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CLDN4\"\u003eCLDN4 function\u003c/a\u003e"}, {"id": "GATA2", "label": "GATA2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=GATA2\"\u003eGATA2 function\u003c/a\u003e"}, {"id": "KLF5", "label": "KLF5", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=KLF5\"\u003eKLF5 function\u003c/a\u003e"}, {"id": "GRHL1", "label": "GRHL1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=GRHL1\"\u003eGRHL1 function\u003c/a\u003e"}, {"id": "DAPK2", "label": "DAPK2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=DAPK2\"\u003eDAPK2 function\u003c/a\u003e"}, {"id": "TP63", "label": "TP63", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TP63\"\u003eTP63 function\u003c/a\u003e"}, {"id": "CD164", "label": "CD164", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CD164\"\u003eCD164 function\u003c/a\u003e"}, {"id": "GADD45G", "label": "GADD45G", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=GADD45G\"\u003eGADD45G function\u003c/a\u003e"}, {"id": "NEUROG1", "label": "NEUROG1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=NEUROG1\"\u003eNEUROG1 function\u003c/a\u003e"}, {"id": "CLPTM1", "label": "CLPTM1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CLPTM1\"\u003eCLPTM1 function\u003c/a\u003e"}, {"id": "DSP", "label": "DSP", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=DSP\"\u003eDSP function\u003c/a\u003e"}, {"id": "TFAP2C", "label": "TFAP2C", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TFAP2C\"\u003eTFAP2C function\u003c/a\u003e"}, {"id": "ARRDC3", "label": "ARRDC3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ARRDC3\"\u003eARRDC3 function\u003c/a\u003e"}, {"id": "TRIM29", "label": "TRIM29", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TRIM29\"\u003eTRIM29 function\u003c/a\u003e"}, {"id": "TSPO", "label": "TSPO", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TSPO\"\u003eTSPO function\u003c/a\u003e"}, {"id": "EPPK1", "label": "EPPK1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=EPPK1\"\u003eEPPK1 function\u003c/a\u003e"}, {"id": "KRT12", "label": "KRT12", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=KRT12\"\u003eKRT12 function\u003c/a\u003e"}, {"id": "EIF1B", "label": "EIF1B", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=EIF1B\"\u003eEIF1B function\u003c/a\u003e"}, {"id": "RREB1", "label": "RREB1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=RREB1\"\u003eRREB1 function\u003c/a\u003e"}, {"id": "CMIP", "label": "CMIP", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CMIP\"\u003eCMIP function\u003c/a\u003e"}, {"id": "ABO", "label": "ABO", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ABO\"\u003eABO function\u003c/a\u003e"}, {"id": "FOXI1", "label": "FOXI1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=FOXI1\"\u003eFOXI1 function\u003c/a\u003e"}, {"id": "TFAP2A", "label": "TFAP2A", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TFAP2A\"\u003eTFAP2A function\u003c/a\u003e"}, {"id": "KLF4", "label": "KLF4", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=KLF4\"\u003eKLF4 function\u003c/a\u003e"}, {"id": "KRT19", "label": "KRT19", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=KRT19\"\u003eKRT19 function\u003c/a\u003e"}, {"id": "GATA3", "label": "GATA3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=GATA3\"\u003eGATA3 function\u003c/a\u003e"}, {"id": "ZEB2", "label": "ZEB2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ZEB2\"\u003eZEB2 function\u003c/a\u003e"}, {"id": "ANXA4", "label": "ANXA4", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ANXA4\"\u003eANXA4 function\u003c/a\u003e"}, {"id": "SPTLC3", "label": "SPTLC3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SPTLC3\"\u003eSPTLC3 function\u003c/a\u003e"}, {"id": "AHNAK", "label": "AHNAK", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=AHNAK\"\u003eAHNAK function\u003c/a\u003e"}, {"id": "LHX5", "label": "LHX5", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=LHX5\"\u003eLHX5 function\u003c/a\u003e"}, {"id": "SP5", "label": "SP5", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SP5\"\u003eSP5 function\u003c/a\u003e"}, {"id": "AXIN2", "label": "AXIN2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=AXIN2\"\u003eAXIN2 function\u003c/a\u003e"}, {"id": "TCF7L1", "label": "TCF7L1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TCF7L1\"\u003eTCF7L1 function\u003c/a\u003e"}, {"id": "PITX1", "label": "PITX1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=PITX1\"\u003ePITX1 function\u003c/a\u003e"}, {"id": "BTG1", "label": "BTG1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=BTG1\"\u003eBTG1 function\u003c/a\u003e"}, {"id": "ACTB", "label": "ACTB", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ACTB\"\u003eACTB function\u003c/a\u003e"}, {"id": "ESRP1", "label": "ESRP1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ESRP1\"\u003eESRP1 function\u003c/a\u003e"}, {"id": "DLL1", "label": "DLL1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=DLL1\"\u003eDLL1 function\u003c/a\u003e"}, {"id": "SALL4", "label": "SALL4", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SALL4\"\u003eSALL4 function\u003c/a\u003e"}, {"id": "CDH1", "label": "CDH1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CDH1\"\u003eCDH1 function\u003c/a\u003e"}, {"id": "SOX11", "label": "SOX11", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SOX11\"\u003eSOX11 function\u003c/a\u003e"}, {"id": "HAS1", "label": "HAS1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=HAS1\"\u003eHAS1 function\u003c/a\u003e"}, {"id": "NR6A1", "label": "NR6A1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=NR6A1\"\u003eNR6A1 function\u003c/a\u003e"}, {"id": "HMGB1", "label": "HMGB1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=HMGB1\"\u003eHMGB1 function\u003c/a\u003e"}, {"id": "SRSF5", "label": "SRSF5", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SRSF5\"\u003eSRSF5 function\u003c/a\u003e"}, {"id": "OTX2", "label": "OTX2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=OTX2\"\u003eOTX2 function\u003c/a\u003e"}, {"id": "OTX1", "label": "OTX1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=OTX1\"\u003eOTX1 function\u003c/a\u003e"}, {"id": "HESX1", "label": "HESX1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=HESX1\"\u003eHESX1 function\u003c/a\u003e"}, {"id": "SRSF7", "label": "SRSF7", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SRSF7\"\u003eSRSF7 function\u003c/a\u003e"}, {"id": "HIST2H2AB", "label": "HIST2H2AB", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=HIST2H2AB\"\u003eHIST2H2AB function\u003c/a\u003e"}, {"id": "CYP2F1", "label": "CYP2F1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYP2F1\"\u003eCYP2F1 function\u003c/a\u003e"}, {"id": "ATP2A2", "label": "ATP2A2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ATP2A2\"\u003eATP2A2 function\u003c/a\u003e"}, {"id": "DLX3", "label": "DLX3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=DLX3\"\u003eDLX3 function\u003c/a\u003e"}, {"id": "SSR3", "label": "SSR3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SSR3\"\u003eSSR3 function\u003c/a\u003e"}, {"id": "LRWD1", "label": "LRWD1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=LRWD1\"\u003eLRWD1 function\u003c/a\u003e"}, {"id": "ANXA2", "label": "ANXA2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ANXA2\"\u003eANXA2 function\u003c/a\u003e"}, {"id": "SEC61A1", "label": "SEC61A1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SEC61A1\"\u003eSEC61A1 function\u003c/a\u003e"}, {"id": "TBX3", "label": "TBX3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TBX3\"\u003eTBX3 function\u003c/a\u003e"}, {"id": "PCDH7", "label": "PCDH7", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=PCDH7\"\u003ePCDH7 function\u003c/a\u003e"}, {"id": "HNRNPDL", "label": "HNRNPDL", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=HNRNPDL\"\u003eHNRNPDL function\u003c/a\u003e"}, {"id": "MKRN1", "label": "MKRN1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=MKRN1\"\u003eMKRN1 function\u003c/a\u003e"}, {"id": "CNDP1", "label": "CNDP1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CNDP1\"\u003eCNDP1 function\u003c/a\u003e"}, {"id": "ZIC3", "label": "ZIC3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ZIC3\"\u003eZIC3 function\u003c/a\u003e"}, {"id": "ZIC1", "label": "ZIC1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ZIC1\"\u003eZIC1 function\u003c/a\u003e"}, {"id": "ACTN1", "label": "ACTN1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ACTN1\"\u003eACTN1 function\u003c/a\u003e"}, {"id": "SNAI2", "label": "SNAI2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SNAI2\"\u003eSNAI2 function\u003c/a\u003e"}, {"id": "MYC", "label": "MYC", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=MYC\"\u003eMYC function\u003c/a\u003e"}, {"id": "SOX9", "label": "SOX9", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SOX9\"\u003eSOX9 function\u003c/a\u003e"}, {"id": "SOX8", "label": "SOX8", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SOX8\"\u003eSOX8 function\u003c/a\u003e"}, {"id": "HES1", "label": "HES1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=HES1\"\u003eHES1 function\u003c/a\u003e"}, {"id": "ZMIZ1", "label": "ZMIZ1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ZMIZ1\"\u003eZMIZ1 function\u003c/a\u003e"}, {"id": "SHISA2", "label": "SHISA2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SHISA2\"\u003eSHISA2 function\u003c/a\u003e"}, {"id": "PHF21B", "label": "PHF21B", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=PHF21B\"\u003ePHF21B function\u003c/a\u003e"}, {"id": "ID2", "label": "ID2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ID2\"\u003eID2 function\u003c/a\u003e"}, {"id": "GFPT1", "label": "GFPT1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=GFPT1\"\u003eGFPT1 function\u003c/a\u003e"}, {"id": "SALL2", "label": "SALL2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SALL2\"\u003eSALL2 function\u003c/a\u003e"}, {"id": "SOX3", "label": "SOX3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SOX3\"\u003eSOX3 function\u003c/a\u003e"}, {"id": "FOXA1", "label": "FOXA1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=FOXA1\"\u003eFOXA1 function\u003c/a\u003e"}, {"id": "KRT15", "label": "KRT15", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=KRT15\"\u003eKRT15 function\u003c/a\u003e"}, {"id": "MKL1", "label": "MKL1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=MKL1\"\u003eMKL1 function\u003c/a\u003e"}, {"id": "BMP7", "label": "BMP7", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=BMP7\"\u003eBMP7 function\u003c/a\u003e"}, {"id": "ZNF706", "label": "ZNF706", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ZNF706\"\u003eZNF706 function\u003c/a\u003e"}, {"id": "CFH", "label": "CFH", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CFH\"\u003eCFH function\u003c/a\u003e"}, {"id": "POU5F1", "label": "POU5F1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=POU5F1\"\u003ePOU5F1 function\u003c/a\u003e"}, {"id": "CALR", "label": "CALR", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CALR\"\u003eCALR function\u003c/a\u003e"}, {"id": "CDKN1B", "label": "CDKN1B", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CDKN1B\"\u003eCDKN1B function\u003c/a\u003e"}, {"id": "INSM1", "label": "INSM1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=INSM1\"\u003eINSM1 function\u003c/a\u003e"}, {"id": "GMNN", "label": "GMNN", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=GMNN\"\u003eGMNN function\u003c/a\u003e"}, {"id": "CLDN6", "label": "CLDN6", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CLDN6\"\u003eCLDN6 function\u003c/a\u003e"}, {"id": "TAGLN2", "label": "TAGLN2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TAGLN2\"\u003eTAGLN2 function\u003c/a\u003e"}, {"id": "FLNA", "label": "FLNA", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=FLNA\"\u003eFLNA function\u003c/a\u003e"}, {"id": "CDX4", "label": "CDX4", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CDX4\"\u003eCDX4 function\u003c/a\u003e"}, {"id": "NGFR", "label": "NGFR", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=NGFR\"\u003eNGFR function\u003c/a\u003e"}, {"id": "RGMA", "label": "RGMA", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=RGMA\"\u003eRGMA function\u003c/a\u003e"}, {"id": "NECTIN1", "label": "NECTIN1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=NECTIN1\"\u003eNECTIN1 function\u003c/a\u003e"}, {"id": "TRIL", "label": "TRIL", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TRIL\"\u003eTRIL function\u003c/a\u003e"}, {"id": "S100A10", "label": "S100A10", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=S100A10\"\u003eS100A10 function\u003c/a\u003e"}, {"id": "SOX17", "label": "SOX17", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SOX17\"\u003eSOX17 function\u003c/a\u003e"}, {"id": "CHP1", "label": "CHP1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CHP1\"\u003eCHP1 function\u003c/a\u003e"}, {"id": "HIC2", "label": "HIC2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=HIC2\"\u003eHIC2 function\u003c/a\u003e"}, {"id": "PTBP1", "label": "PTBP1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=PTBP1\"\u003ePTBP1 function\u003c/a\u003e"}, {"id": "CRABP2", "label": "CRABP2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CRABP2\"\u003eCRABP2 function\u003c/a\u003e"}, {"id": "GJB1", "label": "GJB1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=GJB1\"\u003eGJB1 function\u003c/a\u003e"}, {"id": "S100A11", "label": "S100A11", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=S100A11\"\u003eS100A11 function\u003c/a\u003e"}, {"id": "KRT8", "label": "KRT8", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=KRT8\"\u003eKRT8 function\u003c/a\u003e"}, {"id": "CLCA4", "label": "CLCA4", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CLCA4\"\u003eCLCA4 function\u003c/a\u003e"}, {"id": "TACSTD2", "label": "TACSTD2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TACSTD2\"\u003eTACSTD2 function\u003c/a\u003e"}, {"id": "UGDH", "label": "UGDH", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=UGDH\"\u003eUGDH function\u003c/a\u003e"}, {"id": "DBN1", "label": "DBN1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=DBN1\"\u003eDBN1 function\u003c/a\u003e"}, {"id": "LHX2", "label": "LHX2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=LHX2\"\u003eLHX2 function\u003c/a\u003e"}, {"id": "GLDC", "label": "GLDC", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=GLDC\"\u003eGLDC function\u003c/a\u003e"}, {"id": "PYGB", "label": "PYGB", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=PYGB\"\u003ePYGB function\u003c/a\u003e"}, {"id": "PCNA", "label": "PCNA", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=PCNA\"\u003ePCNA function\u003c/a\u003e"}, {"id": "MYH9", "label": "MYH9", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=MYH9\"\u003eMYH9 function\u003c/a\u003e"}, {"id": "TLL2", "label": "TLL2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TLL2\"\u003eTLL2 function\u003c/a\u003e"}, {"id": "NOTCH1", "label": "NOTCH1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=NOTCH1\"\u003eNOTCH1 function\u003c/a\u003e"}, {"id": "EBF3", "label": "EBF3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=EBF3\"\u003eEBF3 function\u003c/a\u003e"}, {"id": "POU4F3", "label": "POU4F3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=POU4F3\"\u003ePOU4F3 function\u003c/a\u003e"}, {"id": "ISL1", "label": "ISL1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ISL1\"\u003eISL1 function\u003c/a\u003e"}, {"id": "GATA4", "label": "GATA4", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=GATA4\"\u003eGATA4 function\u003c/a\u003e"}, {"id": "CNDP2", "label": "CNDP2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CNDP2\"\u003eCNDP2 function\u003c/a\u003e"}, {"id": "CTBS", "label": "CTBS", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CTBS\"\u003eCTBS function\u003c/a\u003e"}, {"id": "POU3F1", "label": "POU3F1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=POU3F1\"\u003ePOU3F1 function\u003c/a\u003e"}, {"id": "VIM", "label": "VIM", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=VIM\"\u003eVIM function\u003c/a\u003e"}, {"id": "ZIC2", "label": "ZIC2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ZIC2\"\u003eZIC2 function\u003c/a\u003e"}, {"id": "MEIS2", "label": "MEIS2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=MEIS2\"\u003eMEIS2 function\u003c/a\u003e"}, {"id": "CADM1", "label": "CADM1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CADM1\"\u003eCADM1 function\u003c/a\u003e"}, {"id": "KRT18", "label": "KRT18", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=KRT18\"\u003eKRT18 function\u003c/a\u003e"}, {"id": "OLFM4", "label": "OLFM4", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=OLFM4\"\u003eOLFM4 function\u003c/a\u003e"}, {"id": "RRM2", "label": "RRM2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=RRM2\"\u003eRRM2 function\u003c/a\u003e"}, {"id": "PMM1", "label": "PMM1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=PMM1\"\u003ePMM1 function\u003c/a\u003e"}, {"id": "ELAVL3", "label": "ELAVL3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ELAVL3\"\u003eELAVL3 function\u003c/a\u003e"}, {"id": "CAPN8", "label": "CAPN8", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CAPN8\"\u003eCAPN8 function\u003c/a\u003e"}, {"id": "ROR2", "label": "ROR2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ROR2\"\u003eROR2 function\u003c/a\u003e"}, {"id": "C9", "label": "C9", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=C9\"\u003eC9 function\u003c/a\u003e"}, {"id": "PAX8", "label": "PAX8", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=PAX8\"\u003ePAX8 function\u003c/a\u003e"}, {"id": "MEX3B", "label": "MEX3B", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=MEX3B\"\u003eMEX3B function\u003c/a\u003e"}, {"id": "PAX3", "label": "PAX3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=PAX3\"\u003ePAX3 function\u003c/a\u003e"}, {"id": "SERBP1", "label": "SERBP1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SERBP1\"\u003eSERBP1 function\u003c/a\u003e"}, {"id": "RAX", "label": "RAX", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=RAX\"\u003eRAX function\u003c/a\u003e"}, {"id": "PAX6", "label": "PAX6", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=PAX6\"\u003ePAX6 function\u003c/a\u003e"}, {"id": "HNRNPK", "label": "HNRNPK", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=HNRNPK\"\u003eHNRNPK function\u003c/a\u003e"}, {"id": "ZC3H12C", "label": "ZC3H12C", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ZC3H12C\"\u003eZC3H12C function\u003c/a\u003e"}, {"id": "EFNB2", "label": "EFNB2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=EFNB2\"\u003eEFNB2 function\u003c/a\u003e"}, {"id": "C3", "label": "C3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=C3\"\u003eC3 function\u003c/a\u003e"}, {"id": "NKX6-2", "label": "NKX6-2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=NKX6-2\"\u003eNKX6-2 function\u003c/a\u003e"}, {"id": "PTCH1", "label": "PTCH1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=PTCH1\"\u003ePTCH1 function\u003c/a\u003e"}, {"id": "NID2", "label": "NID2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=NID2\"\u003eNID2 function\u003c/a\u003e"}, {"id": "CEP131", "label": "CEP131", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CEP131\"\u003eCEP131 function\u003c/a\u003e"}, {"id": "TET3", "label": "TET3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TET3\"\u003eTET3 function\u003c/a\u003e"}, {"id": "CSGALNACT2", "label": "CSGALNACT2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CSGALNACT2\"\u003eCSGALNACT2 function\u003c/a\u003e"}, {"id": "FOXD3", "label": "FOXD3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=FOXD3\"\u003eFOXD3 function\u003c/a\u003e"}, {"id": "TUBB2B", "label": "TUBB2B", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=TUBB2B\"\u003eTUBB2B function\u003c/a\u003e"}, {"id": "MSI1", "label": "MSI1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=MSI1\"\u003eMSI1 function\u003c/a\u003e"}, {"id": "AKAP12", "label": "AKAP12", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=AKAP12\"\u003eAKAP12 function\u003c/a\u003e"}, {"id": "SRRM4", "label": "SRRM4", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=SRRM4\"\u003eSRRM4 function\u003c/a\u003e"}, {"id": "HS3ST3B1", "label": "HS3ST3B1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=HS3ST3B1\"\u003eHS3ST3B1 function\u003c/a\u003e"}, {"id": "ZBTB16", "label": "ZBTB16", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ZBTB16\"\u003eZBTB16 function\u003c/a\u003e"}, {"id": "NR2F2", "label": "NR2F2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=NR2F2\"\u003eNR2F2 function\u003c/a\u003e"}, {"id": "PRPH", "label": "PRPH", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=PRPH\"\u003ePRPH function\u003c/a\u003e"}, {"id": "PDK1", "label": "PDK1", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=PDK1\"\u003ePDK1 function\u003c/a\u003e"}, {"id": "CELF3", "label": "CELF3", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=CELF3\"\u003eCELF3 function\u003c/a\u003e"}, {"id": "ZC4H2", "label": "ZC4H2", "shape": "dot", "size": 10, "title": "\u003ca href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=ZC4H2\"\u003eZC4H2 function\u003c/a\u003e"}]);
        edges = new vis.DataSet([{"color": "red", "from": "TFCP2L1", "to": "SLC16A3", "weight": 5}, {"color": "red", "from": "TFCP2L1", "to": "ATP6AP2", "weight": 5}, {"color": "red", "from": "TFCP2L1", "to": "BSG", "weight": 5}, {"color": "red", "from": "TFCP2L1", "to": "CLDN4", "weight": 5}, {"color": "red", "from": "TFCP2L1", "to": "GRHL3", "weight": 5}, {"color": "red", "from": "TFCP2L1", "to": "GRHL1", "weight": 5}, {"color": "red", "from": "TFCP2L1", "to": "DAPK2", "weight": 5}, {"color": "red", "from": "TFCP2L1", "to": "CD164", "weight": 5}, {"color": "red", "from": "TFCP2L1", "to": "GADD45G", "weight": 5}, {"color": "red", "from": "TFCP2L1", "to": "CLPTM1", "weight": 5}, {"color": "red", "from": "TFCP2L1", "to": "DSP", "weight": 5}, {"color": "red", "from": "TFCP2L1", "to": "ARRDC3", "weight": 5}, {"color": "red", "from": "SLC16A3", "to": "GRHL3", "weight": 5}, {"color": "red", "from": "GRHL3", "to": "GATA2", "weight": 5}, {"color": "red", "from": "GRHL3", "to": "KLF5", "weight": 5}, {"color": "red", "from": "GATA2", "to": "DAPK2", "weight": 5}, {"color": "red", "from": "GATA2", "to": "DSP", "weight": 5}, {"color": "red", "from": "GATA2", "to": "TRIM29", "weight": 5}, {"color": "red", "from": "GATA2", "to": "TSPO", "weight": 5}, {"color": "red", "from": "GATA2", "to": "EPPK1", "weight": 5}, {"color": "red", "from": "GATA2", "to": "KRT12", "weight": 5}, {"color": "red", "from": "GATA2", "to": "EIF1B", "weight": 5}, {"color": "red", "from": "GATA2", "to": "RREB1", "weight": 5}, {"color": "red", "from": "GATA2", "to": "CMIP", "weight": 5}, {"color": "red", "from": "GATA2", "to": "ABO", "weight": 5}, {"color": "red", "from": "GATA2", "to": "FOXI1", "weight": 5}, {"color": "red", "from": "GATA2", "to": "KRT19", "weight": 5}, {"color": "red", "from": "GATA2", "to": "ANXA4", "weight": 5}, {"color": "red", "from": "GATA2", "to": "SPTLC3", "weight": 5}, {"color": "red", "from": "GATA2", "to": "KLF5", "weight": 5}, {"color": "red", "from": "GATA2", "to": "AHNAK", "weight": 5}, {"color": "red", "from": "GATA2", "to": "KLF4", "weight": 5}, {"color": "red", "from": "GATA2", "to": "ESRP1", "weight": 5}, {"color": "red", "from": "GATA2", "to": "CDH1", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "SOX11", "weight": 5}, {"color": "red", "from": "GATA2", "to": "HAS1", "weight": 5}, {"color": "red", "from": "GATA2", "to": "CYP2F1", "weight": 5}, {"color": "red", "from": "GATA2", "to": "DLX3", "weight": 5}, {"color": "red", "from": "GATA2", "to": "SSR3", "weight": 5}, {"color": "red", "from": "GATA2", "to": "SEC61A1", "weight": 5}, {"color": "red", "from": "GATA2", "to": "TP63", "weight": 5}, {"color": "red", "from": "GATA2", "to": "PCDH7", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "HNRNPDL", "weight": 5}, {"color": "red", "from": "GATA2", "to": "CNDP1", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "ZIC3", "weight": 5}, {"color": "red", "from": "GATA2", "to": "TFAP2C", "weight": 5}, {"color": "red", "from": "GATA2", "to": "TFAP2A", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "ZIC1", "weight": 5}, {"color": "red", "from": "GATA2", "to": "GATA3", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "ZEB2", "weight": 5}, {"color": "red", "from": "GATA2", "to": "ACTN1", "weight": 5}, {"color": "red", "from": "GATA2", "to": "ID2", "weight": 5}, {"color": "red", "from": "GATA2", "to": "GFPT1", "weight": 5}, {"color": "red", "from": "GATA2", "to": "MKL1", "weight": 5}, {"color": "red", "from": "GATA2", "to": "CALR", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "GMNN", "weight": 5}, {"color": "red", "from": "GATA2", "to": "FLNA", "weight": 5}, {"color": "red", "from": "GATA2", "to": "TRIL", "weight": 5}, {"color": "red", "from": "GATA2", "to": "S100A10", "weight": 5}, {"color": "red", "from": "GATA2", "to": "HIC2", "weight": 5}, {"color": "red", "from": "GATA2", "to": "PTBP1", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "CRABP2", "weight": 5}, {"color": "red", "from": "GATA2", "to": "KRT8", "weight": 5}, {"color": "red", "from": "GATA2", "to": "TACSTD2", "weight": 5}, {"color": "red", "from": "GATA2", "to": "UGDH", "weight": 5}, {"color": "red", "from": "GATA2", "to": "MYH9", "weight": 5}, {"color": "red", "from": "GATA2", "to": "TLL2", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "SOX3", "weight": 5}, {"color": "red", "from": "GATA2", "to": "CTBS", "weight": 5}, {"color": "red", "from": "GATA2", "to": "KRT18", "weight": 5}, {"color": "red", "from": "GATA2", "to": "OLFM4", "weight": 5}, {"color": "red", "from": "GATA2", "to": "PMM1", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "ELAVL3", "weight": 5}, {"color": "red", "from": "GATA2", "to": "CAPN8", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "ROR2", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "MEX3B", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "PAX3", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "ZIC2", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "TET3", "weight": 5}, {"color": "blue", "from": "GATA2", "to": "AKAP12", "weight": 5}, {"color": "red", "from": "DAPK2", "to": "TP63", "weight": 5}, {"color": "red", "from": "TP63", "to": "TRIM29", "weight": 5}, {"color": "red", "from": "TP63", "to": "AHNAK", "weight": 5}, {"color": "red", "from": "TP63", "to": "ACTB", "weight": 5}, {"color": "red", "from": "TP63", "to": "CDH1", "weight": 5}, {"color": "red", "from": "TP63", "to": "ANXA2", "weight": 5}, {"color": "red", "from": "TP63", "to": "TBX3", "weight": 5}, {"color": "red", "from": "TP63", "to": "PCDH7", "weight": 5}, {"color": "red", "from": "TP63", "to": "ACTN1", "weight": 5}, {"color": "red", "from": "TP63", "to": "TAGLN2", "weight": 5}, {"color": "red", "from": "TP63", "to": "FLNA", "weight": 5}, {"color": "red", "from": "TP63", "to": "NECTIN1", "weight": 5}, {"color": "red", "from": "TP63", "to": "TRIL", "weight": 5}, {"color": "red", "from": "TP63", "to": "S100A10", "weight": 5}, {"color": "red", "from": "TP63", "to": "S100A11", "weight": 5}, {"color": "red", "from": "TP63", "to": "CLCA4", "weight": 5}, {"color": "red", "from": "TP63", "to": "MYH9", "weight": 5}, {"color": "red", "from": "TP63", "to": "TLL2", "weight": 5}, {"color": "red", "from": "TP63", "to": "OLFM4", "weight": 5}, {"color": "red", "from": "TP63", "to": "PMM1", "weight": 5}, {"color": "red", "from": "GADD45G", "to": "NEUROG1", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "DLL1", "weight": 5}, {"color": "blue", "from": "NEUROG1", "to": "HIST2H2AB", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "PHF21B", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "CDKN1B", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "INSM1", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "DBN1", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "NOTCH1", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "EBF3", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "POU4F3", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "ISL1", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "ELAVL3", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "ZC3H12C", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "TUBB2B", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "SRRM4", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "PRPH", "weight": 5}, {"color": "red", "from": "NEUROG1", "to": "CELF3", "weight": 5}, {"color": "red", "from": "DSP", "to": "TFAP2C", "weight": 5}, {"color": "red", "from": "TFAP2C", "to": "KRT19", "weight": 5}, {"color": "red", "from": "TFAP2C", "to": "AHNAK", "weight": 5}, {"color": "red", "from": "TFAP2C", "to": "ESRP1", "weight": 5}, {"color": "red", "from": "TFAP2C", "to": "TFAP2A", "weight": 5}, {"color": "blue", "from": "TFAP2C", "to": "ZEB2", "weight": 5}, {"color": "blue", "from": "TFAP2C", "to": "SALL2", "weight": 5}, {"color": "blue", "from": "TFAP2C", "to": "SOX3", "weight": 5}, {"color": "red", "from": "TFAP2C", "to": "HIC2", "weight": 5}, {"color": "blue", "from": "TFAP2C", "to": "PCNA", "weight": 5}, {"color": "red", "from": "TFAP2C", "to": "OLFM4", "weight": 5}, {"color": "blue", "from": "TFAP2C", "to": "RRM2", "weight": 5}, {"color": "blue", "from": "TFAP2C", "to": "PDK1", "weight": 5}, {"color": "blue", "from": "TFAP2C", "to": "ZC4H2", "weight": 5}, {"color": "red", "from": "FOXI1", "to": "TFAP2A", "weight": 5}, {"color": "red", "from": "FOXI1", "to": "KLF4", "weight": 5}, {"color": "red", "from": "FOXI1", "to": "HAS1", "weight": 5}, {"color": "red", "from": "FOXI1", "to": "SRSF5", "weight": 5}, {"color": "red", "from": "FOXI1", "to": "SRSF7", "weight": 5}, {"color": "red", "from": "FOXI1", "to": "LRWD1", "weight": 5}, {"color": "red", "from": "FOXI1", "to": "BMP7", "weight": 5}, {"color": "blue", "from": "FOXI1", "to": "NR6A1", "weight": 5}, {"color": "red", "from": "TFAP2A", "to": "KLF4", "weight": 5}, {"color": "blue", "from": "TFAP2A", "to": "ZEB2", "weight": 5}, {"color": "blue", "from": "KLF4", "to": "ZEB2", "weight": 5}, {"color": "red", "from": "KRT19", "to": "GATA3", "weight": 5}, {"color": "blue", "from": "KRT19", "to": "ZEB2", "weight": 5}, {"color": "red", "from": "GATA3", "to": "AHNAK", "weight": 5}, {"color": "red", "from": "GATA3", "to": "ESRP1", "weight": 5}, {"color": "red", "from": "GATA3", "to": "CDH1", "weight": 5}, {"color": "red", "from": "GATA3", "to": "ANXA2", "weight": 5}, {"color": "red", "from": "GATA3", "to": "PCDH7", "weight": 5}, {"color": "red", "from": "GATA3", "to": "ACTN1", "weight": 5}, {"color": "red", "from": "GATA3", "to": "ID2", "weight": 5}, {"color": "red", "from": "GATA3", "to": "TRIL", "weight": 5}, {"color": "red", "from": "GATA3", "to": "S100A10", "weight": 5}, {"color": "red", "from": "GATA3", "to": "UGDH", "weight": 5}, {"color": "red", "from": "GATA3", "to": "TLL2", "weight": 5}, {"color": "red", "from": "GATA3", "to": "CTBS", "weight": 5}, {"color": "red", "from": "GATA3", "to": "OLFM4", "weight": 5}, {"color": "blue", "from": "ZEB2", "to": "ESRP1", "weight": 5}, {"color": "red", "from": "ZEB2", "to": "SALL4", "weight": 5}, {"color": "red", "from": "ZEB2", "to": "SOX11", "weight": 5}, {"color": "blue", "from": "ZEB2", "to": "HAS1", "weight": 5}, {"color": "red", "from": "ZEB2", "to": "ATP2A2", "weight": 5}, {"color": "red", "from": "ZEB2", "to": "ZMIZ1", "weight": 5}, {"color": "blue", "from": "ZEB2", "to": "GFPT1", "weight": 5}, {"color": "red", "from": "ZEB2", "to": "RGMA", "weight": 5}, {"color": "blue", "from": "ZEB2", "to": "OLFM4", "weight": 5}, {"color": "red", "from": "ZEB2", "to": "ELAVL3", "weight": 5}, {"color": "red", "from": "ZEB2", "to": "ROR2", "weight": 5}, {"color": "red", "from": "ZEB2", "to": "CEP131", "weight": 5}, {"color": "red", "from": "ZEB2", "to": "TET3", "weight": 5}, {"color": "red", "from": "ZEB2", "to": "FOXD3", "weight": 5}, {"color": "red", "from": "ZEB2", "to": "AKAP12", "weight": 5}, {"color": "red", "from": "ZEB2", "to": "PDK1", "weight": 5}, {"color": "red", "from": "ZEB2", "to": "ZC4H2", "weight": 5}, {"color": "blue", "from": "AHNAK", "to": "LHX5", "weight": 5}, {"color": "red", "from": "LHX5", "to": "MKRN1", "weight": 5}, {"color": "red", "from": "LHX5", "to": "BMP7", "weight": 5}, {"color": "red", "from": "LHX5", "to": "OTX1", "weight": 5}, {"color": "red", "from": "LHX5", "to": "ZIC1", "weight": 5}, {"color": "red", "from": "SP5", "to": "AXIN2", "weight": 5}, {"color": "red", "from": "SP5", "to": "NGFR", "weight": 5}, {"color": "red", "from": "SP5", "to": "POU5F1", "weight": 5}, {"color": "red", "from": "SP5", "to": "CDX4", "weight": 5}, {"color": "blue", "from": "AXIN2", "to": "TCF7L1", "weight": 5}, {"color": "red", "from": "TCF7L1", "to": "OTX2", "weight": 5}, {"color": "red", "from": "TCF7L1", "to": "OTX1", "weight": 5}, {"color": "red", "from": "TCF7L1", "to": "HESX1", "weight": 5}, {"color": "red", "from": "PITX1", "to": "BTG1", "weight": 5}, {"color": "red", "from": "PITX1", "to": "CLCA4", "weight": 5}, {"color": "blue", "from": "HAS1", "to": "NR6A1", "weight": 5}, {"color": "blue", "from": "HAS1", "to": "HMGB1", "weight": 5}, {"color": "blue", "from": "NR6A1", "to": "SRSF7", "weight": 5}, {"color": "blue", "from": "NR6A1", "to": "LRWD1", "weight": 5}, {"color": "blue", "from": "NR6A1", "to": "MKRN1", "weight": 5}, {"color": "blue", "from": "NR6A1", "to": "ZNF706", "weight": 5}, {"color": "red", "from": "NR6A1", "to": "TAGLN2", "weight": 5}, {"color": "red", "from": "NR6A1", "to": "CADM1", "weight": 5}, {"color": "red", "from": "NR6A1", "to": "SERBP1", "weight": 5}, {"color": "red", "from": "HMGB1", "to": "HNRNPK", "weight": 5}, {"color": "red", "from": "OTX2", "to": "OTX1", "weight": 5}, {"color": "red", "from": "OTX1", "to": "SHISA2", "weight": 5}, {"color": "red", "from": "OTX1", "to": "CRABP2", "weight": 5}, {"color": "red", "from": "OTX1", "to": "HESX1", "weight": 5}, {"color": "blue", "from": "OTX1", "to": "MEIS2", "weight": 5}, {"color": "red", "from": "HESX1", "to": "RAX", "weight": 5}, {"color": "red", "from": "HESX1", "to": "PAX6", "weight": 5}, {"color": "red", "from": "ZIC3", "to": "GMNN", "weight": 5}, {"color": "red", "from": "ZIC3", "to": "ZIC1", "weight": 5}, {"color": "red", "from": "ZIC3", "to": "ZIC2", "weight": 5}, {"color": "red", "from": "ZIC3", "to": "FOXD3", "weight": 5}, {"color": "red", "from": "ZIC1", "to": "GMNN", "weight": 5}, {"color": "red", "from": "ZIC1", "to": "ZIC2", "weight": 5}, {"color": "red", "from": "SNAI2", "to": "MYC", "weight": 5}, {"color": "red", "from": "SNAI2", "to": "HES1", "weight": 5}, {"color": "blue", "from": "SNAI2", "to": "CLDN6", "weight": 5}, {"color": "red", "from": "SNAI2", "to": "SOX9", "weight": 5}, {"color": "red", "from": "SNAI2", "to": "C3", "weight": 5}, {"color": "red", "from": "MYC", "to": "SOX9", "weight": 5}, {"color": "red", "from": "MYC", "to": "SOX8", "weight": 5}, {"color": "red", "from": "SOX9", "to": "HES1", "weight": 5}, {"color": "red", "from": "SOX9", "to": "CFH", "weight": 5}, {"color": "red", "from": "SOX9", "to": "C9", "weight": 5}, {"color": "red", "from": "SOX9", "to": "SOX8", "weight": 5}, {"color": "red", "from": "SOX9", "to": "PAX3", "weight": 5}, {"color": "red", "from": "SOX9", "to": "C3", "weight": 5}, {"color": "red", "from": "SOX9", "to": "NID2", "weight": 5}, {"color": "red", "from": "SOX8", "to": "C9", "weight": 5}, {"color": "red", "from": "SALL2", "to": "TET3", "weight": 5}, {"color": "red", "from": "SALL2", "to": "MSI1", "weight": 5}, {"color": "red", "from": "SALL2", "to": "AKAP12", "weight": 5}, {"color": "red", "from": "SALL2", "to": "PAX6", "weight": 5}, {"color": "red", "from": "FOXA1", "to": "KRT15", "weight": 5}, {"color": "red", "from": "FOXA1", "to": "GJB1", "weight": 5}, {"color": "blue", "from": "CFH", "to": "POU5F1", "weight": 5}, {"color": "blue", "from": "POU5F1", "to": "PAX6", "weight": 5}, {"color": "red", "from": "CDKN1B", "to": "INSM1", "weight": 5}, {"color": "red", "from": "CDX4", "to": "NGFR", "weight": 5}, {"color": "red", "from": "CDX4", "to": "PYGB", "weight": 5}, {"color": "red", "from": "SOX17", "to": "CHP1", "weight": 5}, {"color": "red", "from": "LHX2", "to": "GLDC", "weight": 5}, {"color": "red", "from": "LHX2", "to": "PCNA", "weight": 5}, {"color": "red", "from": "LHX2", "to": "RRM2", "weight": 5}, {"color": "red", "from": "LHX2", "to": "EFNB2", "weight": 5}, {"color": "red", "from": "LHX2", "to": "RAX", "weight": 5}, {"color": "red", "from": "LHX2", "to": "CSGALNACT2", "weight": 5}, {"color": "red", "from": "LHX2", "to": "PAX6", "weight": 5}, {"color": "red", "from": "LHX2", "to": "HS3ST3B1", "weight": 5}, {"color": "red", "from": "LHX2", "to": "ZBTB16", "weight": 5}, {"color": "red", "from": "EBF3", "to": "TUBB2B", "weight": 5}, {"color": "red", "from": "EBF3", "to": "PRPH", "weight": 5}, {"color": "red", "from": "GATA4", "to": "CNDP2", "weight": 5}, {"color": "red", "from": "POU3F1", "to": "VIM", "weight": 5}, {"color": "red", "from": "C9", "to": "PAX8", "weight": 5}, {"color": "red", "from": "RAX", "to": "EFNB2", "weight": 5}, {"color": "red", "from": "RAX", "to": "PAX6", "weight": 5}, {"color": "red", "from": "PAX6", "to": "CSGALNACT2", "weight": 5}, {"color": "red", "from": "PAX6", "to": "ZBTB16", "weight": 5}, {"color": "red", "from": "NKX6-2", "to": "PTCH1", "weight": 5}, {"color": "red", "from": "NKX6-2", "to": "ZC4H2", "weight": 5}, {"color": "red", "from": "ZBTB16", "to": "NR2F2", "weight": 5}]);

        // adding nodes and edges to the graph
        data = {nodes: nodes, edges: edges};

        var options = {
    "configure": {
        "enabled": false
    },
    "edges": {
        "color": {
            "inherit": true
        },
        "smooth": {
            "enabled": false,
            "type": "continuous"
        }
    },
    "interaction": {
        "dragNodes": true,
        "hideEdgesOnDrag": false,
        "hideNodesOnDrag": false
    },
    "physics": {
        "enabled": true,
        "stabilization": {
            "enabled": true,
            "fit": true,
            "iterations": 1000,
            "onlyDynamicEdges": false,
            "updateInterval": 50
        }
    }
};
        
        

        

        network = new vis.Network(container, data, options);
	 
        
        // make a custom popup
        var popup = document.createElement("div");
        popup.className = 'popup';
        popupTimeout = null;
        popup.addEventListener('mouseover', function () {
            console.log(popup)
            if (popupTimeout !== null) {
                clearTimeout(popupTimeout);
                popupTimeout = null;
            }
        });
        popup.addEventListener('mouseout', function () {
            if (popupTimeout === null) {
                hidePopup();
            }
        });
        container.appendChild(popup);


        // use the popup event to show
        network.on("showPopup", function (params) {
            showPopup(params);
        });

        // use the hide event to hide it
        network.on("hidePopup", function (params) {
            hidePopup();
        });


        // hiding the popup through css
        function hidePopup() {
            popupTimeout = setTimeout(function () { popup.style.display = 'none'; }, 500);
        }

        // showing the popup
        function showPopup(nodeId) {
            // get the data from the vis.DataSet
            var nodeData = nodes.get([nodeId]);
            popup.innerHTML = nodeData[0].title;

            // get the position of the node
            var posCanvas = network.getPositions([nodeId])[nodeId];

            // get the bounding box of the node
            var boundingBox = network.getBoundingBox(nodeId);

            //position tooltip:
            posCanvas.x = posCanvas.x + 0.5 * (boundingBox.right - boundingBox.left);

            // convert coordinates to the DOM space
            var posDOM = network.canvasToDOM(posCanvas);

            // Give it an offset
            posDOM.x += 10;
            posDOM.y -= 20;

            // show and place the tooltip.
            popup.style.display = 'block';
            popup.style.top = posDOM.y + 'px';
            popup.style.left = posDOM.x + 'px';
        }
        


        
        network.on("stabilizationProgress", function(params) {
      		document.getElementById('loadingBar').removeAttribute("style");
	        var maxWidth = 496;
	        var minWidth = 20;
	        var widthFactor = params.iterations/params.total;
	        var width = Math.max(minWidth,maxWidth * widthFactor);

	        document.getElementById('bar').style.width = width + 'px';
	        document.getElementById('text').innerHTML = Math.round(widthFactor*100) + '%';
	    });
	    network.once("stabilizationIterationsDone", function() {
	        document.getElementById('text').innerHTML = '100%';
	        document.getElementById('bar').style.width = '496px';
	        document.getElementById('loadingBar').style.opacity = 0;
	        // really clean the dom element
	        setTimeout(function () {document.getElementById('loadingBar').style.display = 'none';}, 500);
	    });
        

        return network;

    }

    drawGraph();

</script>
</body>
</html>
