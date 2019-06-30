module Export
    function htmlExport(data; title="Title")
        html_h = """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<HTML>
<HEAD>
    <META HTTP-EQUIV="CONTENT-TYPE" CONTENT="text/html; charset=utf-8">
    <TITLE>"""*title*"""</TITLE>
    <META NAME="GENERATOR" CONTENT="ClinicalTrialUtilities">
    <META NAME="Operator" CONTENT="PharmCat">
    <STYLE TYPE="text/css">

    @page { margin-left: 1in; margin-right: 1in; margin-top: 0.5in; margin-bottom: 0.5in }
    P { margin-bottom: 0.08in }

    table {
 border-collapse: collapse;
        min-width: 300px;
        width: auto;
 }
     P.cell {
         margin-left: 0.04in;
         margin-right: 0.04in;
         margin-top: 0.04in;
         widows: 0;
         orphans: 0
     }
     P.pbr {
        margin-bottom: 0in;
        line-height: 0.28in;
        widows: 0;
        orphans: 0
 }
     TD.title {
         background-color: #ffffff;
         border: none;
         padding: 0in
     }
     FONT.title {
         font-size: 9pt;
         color: #010205;
         font-family: Arial,serif;
     }
     FONT.cell {
         font-size: 8pt;
         color: #010205;
         font-family: Arial,serif;
     }
     FONT.comment {
         font-size: 7pt;
         color: #010205;
         font-family: Arial,serif;
     }
     TD.cell {
         min-width: 60px;
         max-width: 100px;
         background-color: #ffffff;
         border-top: 1.00pt solid #152935;
         border-bottom: 1.00pt solid #aeaeae;
         border-left: 1.00pt solid #e0e0e0;
         border-right: 1.00pt solid #e0e0e0;
         padding: 0in
     }
     TD.hcell {
         background-color: #cccccc;
         border-top: 1.00pt solid #152935;
         border-bottom: 1.00pt solid #152935;
         border-left: 1.00pt solid #152935;
         border-right: 1.00pt solid #152935;
         padding: 0in
     }
     TD.cell:first-of-type {
         border-left: 1.00pt solid #152935;
     }
     TD.cell:last-of-type {
         border-right: 1.00pt solid #152935;
     }
     TR.cell:last-of-type > TD {
         border-bottom: 2.00pt solid #152935;
     }

    </STYLE>
</HEAD>
<BODY LANG="en-US" DIR="LTR">"""

        html_f = """</BODY>
</HTML>"""
        html_pbr ="""<P LANG="ru-RU" ALIGN=LEFT STYLE="margin-bottom: 0in; line-height: 0.28in; widows: 0; orphans: 0"><BR></P>"""
        rown   = size(data, 1)
        coln   = size(data, 2)
        cnames = names(data)

        io = open("export.html", "w")
        truncate(io, 0)
        println(io, html_h)

        print(io, """
    <TABLE CELLPADDING=0 CELLSPACING=0>
        <TR CLASS=cell>
            <TD COLSPAN="""*string(coln)*""" CLASS=title>
                <P ALIGN=CENTER CLASS=cell>
                <FONT CLASS=title><B>"""*title*"""</B></FONT></P>
            </TD>
        </TR>""")
        print(io, """
    <TR VALIGN=BOTTOM CLASS=cell>""")
        for c = 1:coln
            print(io, """
        <TD CLASS=hcell>
            <P ALIGN=CENTER CLASS=cell>
            <FONT CLASS=cell>"""*string(cnames[c])*"""</FONT></P>
        </TD>""")
        end
        print(io, """
    </TR>""")
        for r = 1:rown
            print(io,"""
        <TR CLASS=cell>""")
            for c = 1:coln
                print(io,"""
            <TD VALIGN=TOP CLASS=cell>
                <P ALIGN=RIGHT CLASS=cell>
                <FONT CLASS=cell><SPAN LANG="ru-RU">"""*string(cellformat(data[r,c]))*"""</SPAN></FONT></P>
            </TD>""")
            end
            print(io,"""
        </TR>""")
        end
        println(io, """
    </TABLE>""")
        print(io, html_f)
        close(io)
    end

    function cellformat(val)
        if isa(val, AbstractFloat)
            return round(val, digits=3)
        else
            return val
        end
    end
end #End module
