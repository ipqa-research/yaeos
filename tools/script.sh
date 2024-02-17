OUTDIR="./tapeout"
mkdir -p "$OUTDIR"
rm "$OUTDIR/*"
rm ./tapeout/*

infile=pr
name=$(basename ${infile})
modflag="YAEOSD"
ADFirstAidKitDIR=../src/adiff/autodiff_api/tapenade #ADFirstAidKit

# Forward
tapenade -tangent -head "ar(arval)/(n, v, t)" \
         -ext ${ADFirstAidKitDIR}/PUSHPOPGeneralLib \
         -tgtmodulename "$modflag" \
         -noisize\
         "${infile}".f90 \
         -O tapeout

# Double forward
tapenade -tangent -head "ar_d(arval_d)/(v, t)" \
         -ext ${ADFirstAidKitDIR}/PUSHPOPGeneralLib  \
         -tgtmodulename "$modflag" \
         -noisize \
         "tapeout/${name}_d.f90" \
         -O tapeout

# Triple forward
tapenade -tangent -head "ar_d_d(arval_d_d)/(v, t)" \
         -ext ${ADFirstAidKitDIR}/PUSHPOPGeneralLib  \
         -tgtmodulename "$modflag" \
         -noisize \
         "tapeout/${name}_d_d.f90" \
         -O tapeout

# Forward-Backward
tapenade -reverse -head "ar_d/(n, nd, v, vd, t, td)" \
         -noisize \
         -ext ${ADFirstAidKitDIR}/PUSHPOPGeneralLib  \
         -adjmodulename "$modflag" \
         "tapeout/${name}_d_d_d.f90" \
         -O tapeout

# Single backward
tapenade -reverse -head "ar(arval)/(n, v, t)" \
         -noisize \
         -ext ${ADFirstAidKitDIR}/PUSHPOPGeneralLib  \
         -adjmodulename "$modflag" \
         "tapeout/${name}_d_d_d_b.f90" \
         -O tapeout


rm tapeout/*msg
rm tapeout/${name}_d.f90
rm tapeout/${name}_d_d.f90
rm tapeout/${name}_d_d_d.f90
rm tapeout/${name}_d_d_d_b.f90
mv tapeout/${name}_d_d_d_b_b.f90 tapeout/${name}_diff.f90

sed -i "s/$modflag//g" tapeout/${name}_diff.f90
sed -i "s/TYPE(UNKNOWNTYPE).*//g" tapeout/${name}_diff.f90
sed -i "s/model%\(.*\) \(=\).*/model%\1 => \1/g" tapeout/${name}_diff.f90
sed -i 's/REAL\*8/REAL(8)/' tapeout/${name}_diff.f90

cp tapeout/${name}_diff.f90 ../example/taperobinson.f90
# findent < tapeout/${name}_diff.f90 > tapeout/${name}_diff.f90

# lfortran fmt -i tapeout/${infile}_diff.f90
# =============================================================================


# Tapenade 3.16 (develop) - 13 Sep 2023 12:36 - Java 17.0.9 Linux
# @@ TAPENADE_HOME=/home/ruther/downs/tapenade_3.16/bin/..
#  Builds a differentiated program.
#  Usage: tapenade [options]* filenames
#   options:
#    -head, -root <proc>     set the differentiation root procedure(s)
#                            See FAQ for refined invocation syntax, e.g.
#                            independent and dependent arguments, multiple heads...
#    -tangent, -d            differentiate in forward/tangent mode (default)
#    -reverse, -b            differentiate in reverse/adjoint mode
#    -vector, -multi         turn on "vector" mode (i.e. multi-directional)
#    -specializeactivity <unit_names or %all%>  Allow for several activity patterns per routine
#    -primal, -p             turn off differentiation. Show pointer destinations
#    -output, -o <file>      put all generated code into a single <file>
#    -splitoutputfiles       split generated code, one file per top unit
#    -outputdirectory, -O <directory>  put all generated files in <directory> (default: .)
#    -I <includePath>        add a new search path for include files
#    -tgtvarname <str>       set extension for tangent variables  (default %d)
#    -tgtfuncname <str>      set extension for tangent procedures (default %_d)
#    -tgtmodulename <str>    set extension for tangent modules and types (default %_diff)
#    -adjvarname <str>       set extension for adjoint variables  (default %b)
#    -adjfuncname <str>      set extension for adjoint procedures (default %_b)
#    -adjmodulename <str>    set extension for adjoint modules and types (default %_diff)
#    -modulename <str>       set extension for tangent&adjoint modules and types (default %_diff)
#    -inputlanguage <lang>   language of  input files (fortran, fortran90,
#                            fortran95, or C)
#    -outputlanguage <lang>  language of output files (fortran, fortran90,
#                            fortran95, or C)
#    -ext <file>             incorporate external library description <file>
#    -nolib                  don't load standard libraries descriptions
#    -i<n>                   count <n> bytes for an integer (default -i4)
#    -r<n>                   count <n> bytes for a real (default -r4)
#    -dr<n>                  count <n> bytes for a double real (default -dr8)
#    -p<n>                   count <n> bytes for a pointer (default -p8)
#    -fixinterface           don't use activity to filter user-given (in)dependent vars
#    -noinclude              inline include files
#    -debugTGT               insert instructions for debugging tangent mode
#    -debugADJ               insert instructions for debugging adjoint mode
#    -tracelevel <n>         set the level of detail of trace milestones
#    -msglevel <n>           set the level of detail of error messages
#    -msginfile              insert error messages in output files
#    -dump <file>            write a dump <file>
#    -html                   display results in a web browser
#    -nooptim <str>          turn off optimization <str> (in {activity, difftypes,
#                            diffarguments, stripprimalmodules, spareinit, splitdiff, 
#                            mergediff, saveonlyused, tbr, snapshot, diffliveness,
#                            deadcontrol, recomputeintermediates,
#                            everyoptim}
#    -version                display Tapenade version information
#  Report bugs to <tapenade@inria.fr>.
