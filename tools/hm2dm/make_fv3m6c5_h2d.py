import os
import os.path

year='2013'
monbeg='07'
daybeg='01'
monend='07'
dayend='05'

ym=year+monbeg
ymd=year+monbeg+daybeg
ftail='.01.'+ymd+'00'

dirsrc="/scratch3/NCEPDEV/stmp2/Denise.Worthen/Bench"+ym+"/gfs."+ymd+"/00"
#print(dirsrc)

mnendn=[]
mnendl=[]
dayfilelist=[]

mnendn=(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
mnendl=(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
#for ii in range(0,len(mnendn)):
#  print('month '+str(ii)+' last day '+str(mnendn[ii]))

#list ice files and create an ncra command to convert hourly files to daily means
#since timestamp goes as 'valid at', need to mess with day number in timestamp
hrlist=['06','12','18','00']
for root, subdirs, files in os.walk(dirsrc):
 for name in files:
  if name.startswith('ice'+year) and name.endswith('06'+ftail+'.nc'):
   year  = str(name[3:7])
   yrmon = str(name[3:9])
   mnend = mnendn
   if(int(year)%4 == 0):mnend=mnendl

   day   = str(name[9:11])
   iday  = int(day)
   mon   = str(name[7:9])
   #print('year = '+year+' month = '+mon+' day = '+day)

   imon  = int(mon)
   imonp = imon-1
   imonn = imon+1
   if(imonp ==  0): imonp = 12
   if(imonn == 13): imonn =  1

   #print(year+'  '+mon+' end of month '+str(mnend[imon-1]))

   filelist=[]
   for hr in hrlist:
    if(hr == '00'): iday=iday+1
    #print(str(iday)+' '+str(imon)+'  '+str(mnend[imon-1]))
    if(iday > mnend[imon-1]): 
      iday = 1
      imon = imonn

    date=str(year).zfill(4)+str(imon).zfill(2)+str(iday).zfill(2)
    #print(date)
    hrfilename='ice'+date+hr+ftail+'.nc'
    dyfilename='ice'+yrmon+str(day)+ftail+'.nc'
    dayfilelist.append(dyfilename)

    if os.path.isfile(str(root)+'/'+hrfilename):
     filelist.append(str(root)+'/'+hrfilename)
     # testing
     #filelist.append(hrfilename)
    #print(filelist)

    if len(filelist) == 4:
     strtowrite='ncra '+' '.join(filelist)+' '+str(root)+'/'+dyfilename
     # testing
     #strtowrite='ncra '+' '.join(filelist)+' '+dyfilename
     print(strtowrite)
