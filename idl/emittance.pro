function Qshift, dQsc, Qs, m
dQm=0.0	
if (m gt 0.0) then dQm=-dQsc/2.0+sqrt((dQsc/2.0)^2+(m*Qs)^2) $
  else dQm=-dQsc/2.0-sqrt((dQsc/2.0)^2+(m*Qs)^2)
return, dQm
end

pro emittance;, data_dir
;target_dir ="/d/bhs01/appel/patric/neon_paret/4.17/"   
;data_dir = "/d/bhs01/appel/patric/neon_paret/4.17/" 
;data_dir="/d/bhs01/appel/patric/emittance/on/4/"
;target_dir="/d/bhs01/appel/patric/emittance/on/4/"     
data_dir="/d/bhs01/appel/patric/proton_56/"
target_dir=data_dir
!P.MULTI=[0,1,2]
!P.CHARSIZE=1.7
!P.CHARTHICK=2.0
!P.THICK=5.0
!Y.THICK=5.0
!X.THICK=5.0
!P.FONT = 0
set_font= "helvetica"     
     
;Eingabeparameter:

druck=1                         ; Ausdruck ja/nein (1/0)

emittance=1
envelope=0
fft_envelope=0
momentum_spread=0
loss=0
phaseadvance=0
offset=0
fft_offset=0
entropy=0
xz_mode=0
fft_xz=0       

fileName='mti.dat'
openr, lun, data_dir+fileName, /get_lun
tem=1.0
pmax=long(0)
while(not eof(lun)) do begin
	readf,lun,tem
	pmax=pmax+1
endwhile     
point_lun,lun,0
data =findgen(pmax*3)
	readf, lun,data
free_lun,lun   
turn2=findgen(pmax)
int2=findgen(pmax)  
loss2=findgen(pmax)
for i=0L, pmax-1, 1 do begin
  turn2[i] = data[3*i]
  int2[i] = data[1+3*i] 
  loss2[i] = data[2+3*i]
endfor

;IDL device commands:

loadct,12

set_plot,"x"
if(druck EQ 1) then begin
  set_plot,"ps"
  device,filename=target_dir+"emittance.eps",/color,/encapsulated,/inches,ysize=6.4 ;9.6
endif 

; read parameter file idl.dat:

openr, 1, data_dir+"idl.dat"
readf,1,numprocs                ; Prozesse
readf,1,e_kin                   ; kinetische Energie/MeV
readf,1,qb                      ; Ladung/qe
readf,1,mb                      ; Masse/mp
readf,1,current                 ; Strom      
readf,1,C                       ; Ringumfang
readf,1,Nelements               ; elements in cell 
readf,1,cell_length             ; Zellenlaenge
readf,1,Np                      ; Teilchen fuer Plot
readf,1,pipe_radius
readf,1,NX
readf,1,NY
close,1

; number of cells in ring:

Ncell=C/cell_length

; Constants:

qe=1.6022e-19
mp=1.6605e-27   
rp=1.55e-18
I0=current
pi=3.14159265359
clight=2.99e8
gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) 
beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0))    
freq0=beta0*clight/C      
N0=current*C/(qe*beta0*clight)   
   
; Oeffnen der Datenfiles:
num_p=long(0)
openr, 1, data_dir+"patric.dat"
tem=fltarr(18)   
while(not EOF(1)) do begin
  readf,1,tem
  num_p=num_p+1
end
close,1

openr, 1, data_dir+"patric.dat"
ppic=fltarr(20,num_p) 
readf,1,ppic
close,1

; Grid:

s=fltarr(num_p)
s(0:num_p-1)=ppic(0,0:num_p-1)

; Emittance:

ex=fltarr(num_p)
ey=fltarr(num_p)
                 

ex(0:num_p-1)=ppic(1,0:num_p-1)
ey(0:num_p-1)=ppic(2,0:num_p-1)   
                   
;number of particles
Int=fltarr(num_p) 
Int_inj=fltarr(num_p)
Int(0:num_p-1)=ppic(10,0:num_p-1) 
NPIC=N0/Int(0)  
Int_inj(0:num_p-1)=ppic(18,0:num_p-1)


; offsets

offsetx=fltarr(num_p)
offsety=fltarr(num_p)
pickupx=fltarr(num_p)
pickupy=fltarr(num_p)

offsetx(0:num_p-1)=ppic(13,0:num_p-1)
offsety(0:num_p-1)=ppic(14,0:num_p-1)
pickupx(0:num_p-1)=ppic(16,0:num_p-1)
pickupy(0:num_p-1)=ppic(17,0:num_p-1)

; mean (RMS) radius:

xrms=fltarr(num_p)
yrms=fltarr(num_p)

xrms(0:num_p-1)=ppic(5,0:num_p-1)
yrms(0:num_p-1)=ppic(6,0:num_p-1)

; btf dipole kicks:

han=hanning(num_p,/double)

dtheta=fltarr(num_p)
dtheta(0:num_p-1)=ppic(15,0:num_p-1)
dtheta_fft=fft(han*dtheta,/double)
maxbtf_fft=max(abs(dtheta_fft)^2)

; FFT of envelopes/offsets:

xrms_fft=fft(xrms)
yrms_fft=fft(yrms)
xmy_fft=fft(hanning(num_p)*(xrms-yrms))
xpy_fft=0.5*fft(hanning(num_p)*(xrms+yrms))
offsetx_fft=fft(han*offsetx,/double)
offsety_fft=fft(han*offsety,/double)
pickupx_fft=fft(han*pickupx,/double)
pickupy_fft=fft(han*pickupy,/double)

dt=(s(1)-s(0))/(beta0*clight)   ;  time step in track
freq=findgen(num_p/2+1) / (num_p*dt) ; resulting frequency grid
tunevec=freq/(beta0*clight/C)

offsetx_fft_tot=0.0
for j=0L, num_p/2 do begin
  offsetx_fft_tot=offsetx_fft_tot+abs(offsetx_fft[j])^2/dt
endfor

print, "Total power:", offsetx_fft_tot

; impedance function

Qr=10.0
freqr=10.0*freq0
Rs=1.0
ei=complex(0.0,1.0)
impt=complexarr(num_p/2+1)
for j=0L, num_p/2 do begin
  impt[j]=(freqr/freq[j])*Rs/(1.0-ei*Qr*(freq[j]/freqr-freqr/freq[j])) ; 
endfor

; other momenta:

xz=fltarr(num_p)
x2y2=fltarr(num_p)

xz(0:num_p-1)=ppic(9,0:num_p-1)
x2y2(0:num_p-1)=ppic(8,0:num_p-1)

xz_fft=fft(han*xz,/double)

; max. radius:

xmax=fltarr(num_p)
ymax=fltarr(num_p)

xmax(0:num_p-1)=ppic(3,0:num_p-1)
ymax(0:num_p-1)=ppic(4,0:num_p-1)

; dp/p

dp=fltarr(num_p)
dp(0:num_p-1)=ppic(7,0:num_p-1)

; loss

nparticles=fltarr(num_p)
nparticles(0:num_p-1)=ppic(10,0:num_p-1)

; phase advance

advancex=fltarr(num_p)
advancey=fltarr(num_p)
advancex(0:num_p-1)=ppic(11,0:num_p-1)
advancey(0:num_p-1)=ppic(12,0:num_p-1)

; plot offset:

IF (offset eq 1) then begin
  maxoffset=100.*max(offsetx(0:num_p-1))*3
  baretune=Ncell*97.3/360.0
  plot, s/(beta0*clight)*1.0e3, 100.*offsetx(0:num_p-1),/xstyle,/ystyle,xtitle="t [ms]", $
        ytitle="offset [cm]", yrange=[-maxoffset,maxoffset]
  oplot, s/(beta0*clight)*1.0e3, 100.*offsety(0:num_p-1), color=150
;oplot,s/(beta0*clight)*1.0e3,0.0002*exp(s/(beta0*clight)*1.0e3/0.0625),color=50,line=2
;plot,s/C,180.0/!PI*2.0*acos((0.5*abs(cos(2.0*!PI*1.0*btftune*s/C)+offsetx(0:num_p-1)/maxoffset))),/xstyle,/ystyle,xtitle="turns",ytitle="offset in x",yrange=[0,180]
;plot, s/(beta0*clight)*1.0e3, dtheta(0:num_p-1),/xstyle,/ystyle,xtitle="t [ms]",ytitle="offset in x"
endif

if( entropy eq 1) then begin
  plot, s/(beta0*clight)*1.0e3,dtheta[0:num_p-1],/xstyle,/ystyle,xtitle="t [ms]",title="Entropy"
endif

if( fft_offset eq 1) then begin
  baretune=Ncell*97.3/360.0
  syntune=8.2e3/freq0  
  dQsc=0.009/2.0*0.0
  maxoffset=1.2*1.0e7*max(abs(offsetx_fft(0:num_p/2))^2)
  maxtune=freq[!C]/(beta0*clight/C)
  minoffset=0.0*maxoffset
  plot,tunevec(0:num_p/2),1.0e7*abs(offsetx_fft(0:num_p/2))^2,$ /ylog,$
       /xstyle,/ystyle,title="",xtitle="f/f!L0!N",ytitle="P [arb.units]", $
       yrange=[minoffset,maxoffset],xrange=[baretune-0.015,baretune+0.015]
  print, baretune
  print, maxtune
  plotS, [baretune,baretune],[minoffset,maxoffset],line=1
  plotS, [baretune-syntune,baretune-syntune],[minoffset,maxoffset],line=2,color=200
  plotS, [baretune+syntune,baretune+syntune],[minoffset,maxoffset],line=2,color=200
  plotS, [baretune+Qshift(dQsc,syntune,1),baretune+Qshift(dQsc,syntune,1)], $
         [minoffset,1.0*maxoffset],line=2,color=100
  plotS, [baretune+Qshift(dQsc,syntune,-1),baretune+Qshift(dQsc,syntune,-1)], $
         [minoffset,1.0*maxoffset],line=2,color=100
endif

;plot emittance:
                                      

if( emittance eq 1) then begin   
	;plot, s/cell_length, Int*NPIC,  xtitle="turns", ytitle="Intensity", title=data_dir 
   
	;plot, s/cell_length, loss3, xtitle="turns", ytitle="Loss [%]"     
	plot, s/cell_length, (1-Int/Int_inj)*100  
	                  	
	
	plot,s/cell_length,ex*4/1e-6, xtitle="turns",ytitle=textoidl("\epsilon [mm mrad]"), /nodata
    oplot,s/cell_length,ey*4/1e-6,line=2 
    oplot,s/cell_length,ex*4/1e-6     
	legend, [textoidl('\epsilon_x'), textoidl('\epsilon_y')],line=[0,1], box=0, charsize=1.,/right, /bottom                   
   	;ex=4*ex
	;ey=4*ey   
	;dqx=Int2*NPIC*rp*1.4/(pi*beta0^2*gamma0^3)*1/(ex+sqrt(ex*ey))       
	;dqy=Int2*NPIC*rp*1.4/(pi*beta0^2*gamma0^3)*1/(ey+sqrt(ex*ey))   
	
	;plot, s/cell_length, dqy, xtitle="turns", ytitle=textoidl("|\Delta Q|"), yrange=[0,0.4], /ystyle, /nodata           
	;oplot, s/cell_length, dqx 
	;oplot, s/cell_length, dqy ,line=2
	;hline, 0.5*0.36  
	;legend, [textoidl('\Delta Q_x'), textoidl('\Delta Q_y')], line=[0,1], box=0, charsize=1., /right
	
endif

;plot envelopes:

IF(envelope eq 1) then begin
  plot, s/cell_length, 120.0*xrms, /xstyle, /ystyle,title="",xtitle="cells", $
       ytitle="x and y [cm]",/nodata
  oplot,s/cell_length,100.0*xrms,color=50
  oplot,s/cell_length,100.0*yrms,color=100
;oplot,s/C,100.0*xrms[0]*exp(s/(40.0*20.0*C)),line=2
endif

;plot ffts:

IF (fft_envelope eq 1) then begin
  baretune=Ncell*97.3/360.0
  plot,tunevec(0:num_p/2),abs(xmy_fft(1:num_p/2))^2,/ylog, $
       /xstyle,/ystyle,title="",xtitle="tune",ytitle="P",xrange=[baretune-1.0,baretune+10.0]
  oplot,tunevec(0:num_p/2),abs(xpy_fft(1:num_p/2))^2
endif

;plot momentum spread:

if (momentum_spread eq 1) then begin
  plot,s/(beta0*clight)*1.0e3,sqrt(5.0)*1.0e2*smooth(dp,5), $
       /xstyle,/ystyle,title="",xtitle="t [ms]",ytitle="dp/p [%]"
endif

;plot losses:

;IF (loss eq 1) then begin
;  plot,s/C,nparticles,/xstyle,/ystyle,title="",xtitle="turns",ytitle="N"
;endif

IF (phaseadvance eq 1) then begin
  plot,s/C,Ncell*advancex/(2.0*pi),$
       /xstyle,/ystyle,title="",xtitle="turns",ytitle="sigma"
  oplot,s/C,Ncell*advancex/(2.0*pi)*advancey
endif

IF (xz_mode eq 1) then begin
  maxxz=max(xz(0:num_p-1))
  plot, s/(beta0*clight)*1.0e3, xz(0:num_p-1),/xstyle,/ystyle,xtitle="t [ms]",ytitle="xz [m*m]"
endif

IF (fft_xz eq 1) then begin
  baretune=Ncell*97.3/360.0
  syntune=8.2e3/freq0  
  dQsc=0.0055*1.0

  maxoffset=2.0*max(abs(xz_fft(0:num_p/2))^2)
  maxtune=freq[!C]/(beta0*clight/C)
  minoffset=1.0e-3*maxoffset

  plot,tunevec(0:num_p/2),abs(xz_fft(0:num_p/2))^2,/ylog,$
       /xstyle,/ystyle,title="",xtitle="f/f!L0!N",ytitle="P [arb.units]", $
       yrange=[minoffset,maxoffset],xrange=[baretune-0.01,baretune+0.01]
  print, baretune
  print, maxtune
  plotS, [baretune,baretune],[minoffset,maxoffset],line=1
  plotS, [baretune-syntune,baretune-syntune],[minoffset,maxoffset],line=2,color=200
  plotS, [baretune+syntune,baretune+syntune],[minoffset,maxoffset],line=2,color=200
  plotS, [baretune+Qshift(dQsc,syntune,1),baretune+Qshift(dQsc,syntune,1)], $
         [minoffset,1.0*maxoffset],line=2,color=100
  plotS, [baretune+Qshift(dQsc,syntune,-1),baretune+Qshift(dQsc,syntune,-1)], $
         [minoffset,1.0*maxoffset],line=2,color=100
endif

;plot,tunevec(0:num_p/2),real_part(impt),/xstyle,xrange=[0.0,20.0],yrange=[-0.6,1.2]
;oplot,tunevec(0:num_p/2),imaginary(impt)

if(druck EQ 1) then begin
  device,/close
  spawn,'gv '+target_dir+'emittance.eps'
  set_plot,"x"
endif

end
