OUTDIR="./tapeout"
mkdir -p "$OUTDIR"
rm "$OUTDIR/*"
rm ./tapeout/*

infile="$1"

# Remove the extension (file.f90 => file)
infile="${infile%.*}"

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



# Remove the intermediate files
rm tapeout/*msg
rm tapeout/${name}_d.f90
rm tapeout/${name}_d_d.f90
rm tapeout/${name}_d_d_d.f90
rm tapeout/${name}_d_d_d_b.f90
mv tapeout/${name}_d_d_d_b_b.f90 tapeout/${name}_diff.f90


# Remove the flags
sed -i "s/$modflag//g" tapeout/${name}_diff.f90

# Change the way reals are defined
sed -i 's/REAL\*8/REAL(pr)/' tapeout/${name}_diff.f90
sed -i 's/_8/_pr/' tapeout/${name}_diff.f90

# Remove the implicit R constant that tapenade can not find
sed -i "s/REAL :: r//g" tapeout/${name}_diff.f90

# make the types into classes for polymorphism
sed -i 's/TYPE(\(.*\)),/class(\1),/' tapeout/${name}_diff.f90

# remove all the external definition of functions. we use interfaces now
sed -i '/EXTERNAL.*/d' tapeout/${name}_diff.f90

# Options to format document
# fprettify -i 3 --case 1 1 1 1
# findent < tapeout/${name}_diff.f90 > tapeout/${name}_diff.f90
# lfortran fmt -i tapeout/${infile}_diff.f90
# =============================================================================