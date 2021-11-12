include("file_access.jl")
using Latexify
using Statistics
##
tablefolder = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/Julia results/"

function write_table(io, ab_name)
    (t1,t2) = make_rows(ab_name)
        write(io, t1)
        write(io, t2)
end

function write_short_table(io, ab_name)
    (t1,t2) = make_rows(ab_name; short = true)
        write(io, t1)
        write(io, t2)
end

function make_rows(ab_name; short = false)
    start = [];
    n = []
    thetas = get_start_theta("all")
    h5open(snpanalysis, "r") do fid
        push!(n,length(fid[ab_name]))
        for site in fid[ab_name]
            rs = read(site["bayes_avg1.0"])
            mut = read(site["mut"])  # [backmut, forwardmut]
            wt_aa = (*)("ZZ",read(site["wt_mut_aa/1"])...,"ZZ")
            mut_aa = (*)("ZZ",read(site["wt_mut_aa/2"])...,"ZZ")
            sigs = vec([ii*jj for (ii,jj) in Iterators.product(rs, thetas)])
            if short == false
                newrow = vcat( read(site["loc"]), wt_aa , mut_aa ,mut..., quantile(rs,(.1,.5,.9))..., quantile(sigs,(.1,.5,.9))...)
            else
                newrow = vcat( read(site["loc"]), wt_aa , mut_aa , mut..., quantile(sigs,(.1,.5,.9))...)
            end
            push!(start, newrow)
        end
    end
    t1= "\\multirow{$(n[1])}{*}{$ab_name}"
    t2=latexify(hcat(start...),
        fmt=FancyNumberFormatter(2), transpose=true, env=:tabular, side = repeat(["-"], n[1]))
    return (t1,t2)
end

##
open("$(tablefolder)shorttab.txt", "w+") do io
    write(io, "     &      &      &               &  \\multicolumn{3}{c}{\$\\sigma\$} \\\\ ")
    write(io, "bNAb & site & \$ \\mu \$ & \$ \\mu^\\dagger \$ &  0.1 & 0.5 & 0.9 \\\\ ")
    for ab in ["10-1074","3BNC117"]
        write_short_table(io, ab)
    end
end
##
open("$(tablefolder)longtab.txt", "w+") do io
    write(io, "     &      & & &     &               & \\multicolumn{3}{c}{\$ \\sigma/\\theta \$} & \\multicolumn{3}{c}{\$\\sigma\$} \\\\ ")
    write(io, "bNAb & site & wt AA & mut AA & \$ \\mu \$ & \$ \\mu^\\dagger \$  & 0.1 & 0.5 & 0.9                         & 0.1 & 0.5 & 0.9 \\\\ ")
    for ab in ablist
        write_table(io, ab)
    end
end

##
function format_tab(filename)
    s = open(filename) do io
        read(io, String)
    end
    s = replace(s, "\$-\$" => "")
    s = replace(s, "\$ZZ" => "")
    s = replace(s, "ZZ\$" => "")
    s = replace(s, "\\begin{tabular}{" => "")
    s = replace(s, r"cc+}" => "") # regex
    s = replace(s, "\\end{tabular}" => "")
    open(filename, "w+") do io
        write(io, s)
    end
end
##
format_tab("$(tablefolder)shorttab.txt")
format_tab("$(tablefolder)longtab.txt")