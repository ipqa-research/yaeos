OUTDIR="./tapeout"
mkdir -p "$OUTDIR"
rm "$OUTDIR/*"
rm ./tapeout/*

infile="$1"

name=$(basename ${infile})
modflag="YAEOSD"
ADFirstAidKitDIR=../../../src/adiff/autodiff_api/tapenade #ADFirstAidKit

# Forward
tapenade -tangent -head "excess_gibbs(ge)/(n, t)" \
         -noisize\
         -ext ${ADFirstAidKitDIR}/PUSHPOPGeneralLib \
         -tgtmodulename "$modflag" \
         "${infile}".f90 \
         -O tapeout

# Double forward
tapenade -tangent -head "excess_gibbs_d(ge_d)/(t)" \
         -ext ${ADFirstAidKitDIR}/PUSHPOPGeneralLib  \
         -tgtmodulename "$modflag" \
         -noisize \
         "tapeout/${name}_d.f90" \
         -O tapeout

# Triple forward
tapenade -tangent -head "excess_gibbs_d_d(ge_d_d)/(t)" \
         -noisize \
         -ext ${ADFirstAidKitDIR}/PUSHPOPGeneralLib  \
         -tgtmodulename "$modflag" \
         "tapeout/${name}_d_d.f90" \
         -O tapeout

# Forward-Backward
tapenade -reverse -head "excess_gibbs_d/(n, nd, t, td)" \
         -noisize \
         -ext ${ADFirstAidKitDIR}/PUSHPOPGeneralLib  \
         -adjmodulename "$modflag" \
         "tapeout/${name}_d_d_d.f90" \
         -O tapeout

# Single backward
tapenade -reverse -head "excess_gibbs(ge)/(n, t)" \
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
sed -i "s/TYPE(NRTL) :: model/CLASS(NRTL) :: model/" tapeout/${name}_diff.f90

# sed -i "s/model%\(.*\) \(=\) \(.*\)/model%\1 => \3/g" tapeout/${name}_diff.f90
sed -i 's/REAL\*8/REAL(pr)/' tapeout/${name}_diff.f90

# remove all the external definition of functions. we use interfaces now
sed -i '/EXTERNAL.*/d' tapeout/${name}_diff.f90