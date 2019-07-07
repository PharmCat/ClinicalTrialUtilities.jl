module Export
    function htmlExport(data; sort=NaN, rspan=:all, title="Title", dict = NaN, file=NaN)
        rowlist = Array{String,1}(undef, 0)
        cnames  = names(data)
        if isa(sort, Array)
            sort = Symbol.(sort)
            filter!(x->x in cnames, sort)
            if length(sort) == 0 sort = [cnames[1]] end
        else
            sort = [cnames[1]]
        end
        if isa(rspan, Array)
            rspan = Symbol.(rspan)
            filter!(x->x in sort, rspan)
            if length(sort) == 0 rspan = NaN end
        else
            rspan = sort
        end
        out = ""
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
        width: 100%;
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
         border-top: 1.00pt solid #aeaeae;
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
     TD.midcell {
        min-width: 60px;
        max-width: 100px;
        background-color: #ffffff;
        border-top: 1.00pt solid #aeaeae;
        border-bottom: 1.00pt solid #aeaeae;
        border-left: 1.00pt solid #e0e0e0;
        border-right: 1.00pt solid #e0e0e0;
        padding: 0in
     }
     TD.foot {
        border-top: 2.00pt solid #152935;
        border-bottom: 2.00pt solid #152935;
     }

    </STYLE>
</HEAD>
<BODY LANG="en-US" DIR="LTR" leftmargin="40" rightmargin="20" topmargin="20" bottommargin="20">"""

        html_f = """</BODY>
</HTML>"""
        html_pbr ="""<P LANG="ru-RU" ALIGN=LEFT STYLE="margin-bottom: 0in; line-height: 0.28in; widows: 0; orphans: 0"><BR></P>"""
        rown        = size(data, 1)
        coln        = size(data, 2)
        tablematrix = zeros(Int, rown, coln)

        out *= html_h

        out *= """
    <TABLE CELLPADDING=0 CELLSPACING=0>
        <THEAD>
        <TR CLASS=cell>
            <TD COLSPAN="""*string(coln)*""" CLASS=title>
                <P ALIGN=CENTER CLASS=cell>
                <FONT CLASS=title><B>"""*title*"""</B></FONT></P>
            </TD>
        </TR>"""

        out *= """
    <TR VALIGN=BOTTOM CLASS=cell>"""

        for c = 1:coln
            out *= """
        <TD CLASS=hcell>
            <P ALIGN=CENTER CLASS=cell>
            <FONT CLASS=cell>"""*string(cnames[c])*"""</FONT></P>
        </TD>"""
        end

        out *= """
    </TR>
    </THEAD>

    <TBODY>"""

        sort!(data, tuple(sort))
        tablematrix .= 1
        for c = 1:coln
            s = true
            while s
                s = false
                for r = 2:rown
                    if tablematrix[r,c] !=0 && data[r,c] == data[r-1,c]
                        tablematrix[r,c] -= 1;
                        tablematrix[r-1,c] += 1;
                        s = true;
                    end
                end
            end
        end
        #print(tablematrix)

        for r = 1:rown
            rowstr = ""
            for c = 1:coln
                if tablematrix[r,c] > 0 || !any(x -> x == cnames[c], rspan)
                    rowstr *= """
            <TD ROWSPAN="""*(any(x -> x == cnames[c], rspan) ? string(tablematrix[r,c]) : "1")*""" VALIGN=TOP CLASS=\""""*((c > 1 && c < coln) ? "midcell" : "cell")*"""\">
                <P ALIGN=RIGHT CLASS=cell>
                <FONT CLASS=cell><SPAN LANG="ru-RU">"""*string(cellformat(data[r,c]))*"""</SPAN></FONT></P>
            </TD>"""
                end
            end
            push!(rowlist, rowstr)
        end

        for r in rowlist
            out *="""
        <TR CLASS=cell>"""*r*"""
        </TR>"""
        end

        out *= """
    </TBODY>
    <TFOOT><TR><TD COLSPAN="""*string(coln)*""" class=foot></TD><TR></TFOOT>
    </TABLE>"""
        out *= html_f

        if isa(file, String)
            try
                io = open(file, "w")
                truncate(io, 0)
                print(io, out)
                close(io)
                return true
            catch
                return false
            end
        else
            return out
        end
    end

    function cellformat(val)
        if isa(val, AbstractFloat)
            return round(val, digits=3)
        else
            return val
        end
    end
end #End module
