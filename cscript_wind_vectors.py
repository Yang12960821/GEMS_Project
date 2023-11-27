from climaf.api import *

test_dict = dict(project='IGCM_OUT',
             login='p24balk',
             model='IPSLCM6',
             experiment='piControl',
             frequency='monthly',
             simulation='CM61-LR-5perc-pi-AER-01',
             period='2000',
             status='PROD',
             )
ref_dict = dict(project='IGCM_OUT',
             login='p24balk',
             model='IPSLCM6',
             experiment='historical',
             frequency='monthly',
             simulation='CM61-LR-NoDust-hist-AER-01',
             period='2000',
             status='PROD',
             )
clog('debug')
ref_vitu = ds(variable='vitu', **ref_dict).explore('resolve')
summary(ref_vitu)
ref_vitv = ds(variable='vitv', **ref_dict).explore('resolve')
summary(ref_vitv)
ref_pres = ds(variable='pres', **ref_dict).explore('resolve')
summary(ref_pres)

test_vitu = ds(variable='vitu', **test_dict).explore('resolve')
test_vitv = ds(variable='vitv', **test_dict).explore('resolve')
test_pres = ds(variable='pres', **test_dict).explore('resolve')

cscript('yves_script','ferret -unmapped -script /home/jservon/C-ESM-EP/pre_v2.0/share/scientific_packages/dust_diags/script_vectors_SurfWinds.jnl ${in_1} ${in_2} ${in_3} ${in_4} ${in_5} ${out} ${title}',
        format='png')
cscript('yves_script2','ferret -unmapped -script /home/jservon/C-ESM-EP/pre_v2.0/share/scientific_packages/dust_diags/script_vectors_SurfWinds.jnl ${in_1} ${in_2} ${in_3} ${in_4} ${in_5} ${title} ; Fprint -o SurfWindsChange.ps metafile.plt ; ps2png SurfWindsChange.ps ; convert SurfWindsChange.png -rotate -90 ${out} ; rm -f SurfWindsChange.ps metafile.plt SurfWindsChange.png ',
        format='png')
cscript('yves_script3','ferret -unmapped -script /home/jservon/C-ESM-EP/pre_v2.0/share/scientific_packages/dust_diags/script_vectors_SurfWinds.jnl ${in_1} ${in_2} ${in_3} ${in_4} ${in_5} ${title} ; Fprint -o SurfWindsChange.ps metafile.plt ; ps2png SurfWindsChange.ps ; cp SurfWindsChange.png ${out} ; rm -f SurfWindsChange.ps metafile.plt SurfWindsChange.png ',
        format='png')
#sp Fprint -o $TMPDIR/SurfWindsChange.ps metafile.plt
# ps2png $TMPDIR/SurfWindsChange.ps
# rm $TMPDIR/SurfWindsChange.ps
# mv $TMPDIR/SurfWindsChange.png "$6"

cfile()
csync(True)

clog('critical')
test = yves_script2(ref_vitu, ref_vitv, ref_pres, test_vitu, test_vitv,
                   title='MySimtest')
cdrop(test)
iplot(test)

!cat last.out
cdrop(test)
cfile(test)
!cat last.out
!ncdump -h /data/jservon/climafcache/b09b1/0cfe9/1ceaa/e5a60/f18ca/b3c20/30ed0/315ea/cbbf5/eb026/d3c6e/a.nc
csync(True)
test = minus(time_average(ds('IGCM_OUT%CM61-LR-5perc-pi-AER-01%tos%1985-2014%global%/ccc/store/cont003/thredds%p24balk%IPSLCM6%PROD%piControl%OCE%Analyse%*%monthly%last_30Y%*')),time_average(ds('IGCM_OUT%CM61-LR-NoDust-hist-AER-01%tos%1985-2014%global%/ccc/store/cont003/thredds%p24balk%IPSLCM6%PROD%historical%OCE%Analyse%*%monthly%last_30Y%*')))
cfile(test)
!cat last.out

implot(test, mpCenterLonF=0, contours=1)

calias('IGCM_OUT','ovap' , filenameVar='histmth')
calias('IGCM_OUT','temp' , filenameVar='histmth')
calias('IGCM_OUT','zfull', filenameVar='histmth')
calias('IGCM_OUT','zhalf', filenameVar='histmth')
#/ccc/store/cont003/gen2201/p24balk/IGCM_OUT/LMDZORINCA/PROD/AER/LOI6.5perc/ATM/Output/MO/LOI6.5perc_19800101_19801231_1M_histmth.nc

#/ccc/store/cont003/thredds/p24balk/LMDZORINCA/PROD/AER/LOI6.5perc/ATM/Output/MO


ref_dict = dict(project='IGCM_OUT',
                root='/ccc/store/cont003/thredds',
             login='p24balk',
             model='LMDZORINCA',
             experiment='AER',
             frequency='monthly',
             simulation='LOI6.5perc',
             OUT='Output',
             period='1980',
             status='PROD',
             )

ref_ovap  = ds(variable='ovap', **ref_dict).explore('resolve')
ref_temp  = ds(variable='temp', **ref_dict).explore('resolve')
ref_zfull = ds(variable='zfull', **ref_dict).explore('resolve')
ref_zhalf = ds(variable='zhalf', **ref_dict).explore('resolve')
test = compute_MSE(ref_ovap, ref_temp, ref_zfull, ref_zhalf)
cdrop(test)
cfile(test)
test_dict

test_ovap  = ds(variable='ovap',  OUT='Output', **test_dict).explore('resolve')
test_temp  = ds(variable='temp',  OUT='Output', **test_dict).explore('resolve')
test_zfull = ds(variable='zfull', OUT='Output', **test_dict).explore('resolve')
test_zhalf = ds(variable='zhalf', OUT='Output', **test_dict).explore('resolve')

MSE_test = compute_MSE(test_ovap, test_temp, test_zfull, test_zhalf)
clog('debug')
cdrop(MSE_test)
cfile(MSE_test)
ncdump(MSE_test)

!rm -f /data/jservon/climafcache/07e89/f413c/e4c51/25951/ace27/a9a25/3ccc3/78d71/2ffca/91337/88fad/1_61971.nc
clog('critical')
plotdiff = plot(minus(MSE_test, test), min=-20000, max=20000, delta=2000, color='BlWhRe',)
iplot(plotdiff)
!cat last.out

ncdump(time_average(ccdo(test, operator='sellevidx,79')))
implot(time_average(ccdo(test, operator='sellevidx,79')),
       color='MPL_Reds', 
       min=200000, max=380000, delta=20000,
       #colors='20 22 24 26 28 30 32 34 36 38',
       #colors='200000 220000 240000 260000 280000 300000 320000 340000 360000 380000',
       contours=1,
       title='My Model',
       gsnLeftString='1980-2005',
       gsnRightString='MSE (J.kg-1)'
      )

#200000 et 380000 (on peut le faire par incrément de 20000

!cat last.out
implot(time_average(test))
!ls /ccc/store/cont003/thredds/p24balk/IPSLCM6/PROD/historical/CM61-LR-NoDust-hist-AER-01/*/Analyse/SE
!ls /ccc/store/cont003/thredds/p24balk/IPSLCM6/PROD/historical/CM61-LR-NoDust-hist-AER-01/*/Analyse/SE*/CM61-LR-NoDust-hist-AER-01_SE_915_1924_1M_histmth.nc


