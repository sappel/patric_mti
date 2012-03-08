pro follow
data_dir ='/d/bhs01/appel/patric/mtj2/'
outpath = '/u/sappel/tmp/'  
   




       
!P.CHARSIZE=1.6*1.5
!P.CHARTHICK=2.0*1.5
!P.THICK=3.5
!Y.THICK=6.5*1.5
!X.THICK=6.5*1.5

pict_max=101
loadct,39                       ; Color: 39 (Rainbow), 3 (red temperature)
ds=0.4                          ; determines xsmax and ysmax
smoothing=2                     ; contour smoothing
clevels=40                      ; contour levels     

;plot ps files
tvlct, 255, 0, 0, 10 ;rot
tvlct, 0, 0, 255, 11 ;blau
tvlct, 0, 255, 0, 12 ;gruen

;set_plot,"x"
;device,decomposed=0
xs_ps = 8.8
ys_ps = 5.5
x_size=640
y_size=512

set_plot,"ps"
!p.font=0                       ; for font selection
;  device, /helvetica, xsize=x_size, ysize=y_size
;destPath = outpath
cs = 1.5                        ; character size

; read parameter file idl.dat:
openr,1, data_dir+"idl.dat"
readf,1, numprocs                ; Prozesse
readf,1, e_kin                   ; kinetische Energie/MeV
readf,1, qb                      ; Ladung/qe
readf,1, mb                      ; Masse/mp
readf,1, current                 ; Strom      
readf,1, C                       ; Ringumfang
readf,1, Nelements               ; elements in cell 
readf,1, cell_length             ; Zellenlaenge
readf,1, Np                      ; Teilchen fuer Plot
readf,1, pipe_radius
readf,1, NX
readf,1, NY
readf,1, NZ
readf,1, print_cell     
readf,1, Q_hor
readf,1, Q_vec
readf,1, betax
readf,1, alpx
readf,1, x0
readf,1, x0prime   
readf,1, amp0    
readf,1, ampp0
readf,1, delAmp  
readf,1, cells 
readf,1, max_inj
close,1     

;amp0=0
print, 'amp0:', amp0/1e-2 
print,  'ampp0:', ampp0/1e-3 
print, 'x0:', x0/1e-2
print, 'x0prime:', x0prime/1e-3      
           
qe=1.6022e-19
mp=1.6605e-27
I0=current
pi=3.14159265359
clight=2.99e8
gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) 
beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0))   

xsmax=1.0e3*4.0*pipe_radius/(ds*(NX-1.0))
xs=fltarr(NX)
xs(0)=-xsmax
for i=1, NX-1 do xs(i)=xs(i-1)+2.0*xsmax/(NX-1.0)

pics=fltarr(pict_max, 8, Np*numprocs)

x = findgen(NX)*2*pipe_radius/(NX-1)-pipe_radius

; Oeffnen der PIC Datenfiles:

num_p=0L
openr, 1, data_dir+"pics_0.dat"
tem=fltarr(8,Np)   
while(not EOF(1)) do begin
  readu,1,tem
  pics(num_p,0:7,0:Np-1)=tem(0:7,0:Np-1)
  num_p=num_p+1
end
close,1

for i=1,numprocs-1 do begin
  openr, 1, data_dir+"pics_"+string(strtrim(i,1))+".dat"   
  for j=0, num_p-1 do begin
    readu,1,tem
    pics(j,0:7,i*Np:(i+1)*Np-1)=tem(0:7,0:Np-1)
  endfor
  close,1 
endfor

rho_x=fltarr(num_p,NX)
rho_xy=fltarr(num_p,NX,NY)


IF 0 THEN BEGIN ; plot xprofile?
for j=0,num_p-1 do BEGIN
  !p.title=strcompress("time = " + $
                       string(format='(F7.2)', 1e6*j*cell_length*print_cell/(beta0*clight)) $
                       + " us")

  if (j lt 9) then begin
    device,filename=strcompress(outpath+"rhox00" + string(j+1) + ".eps", /remove_all), $
           /color, /encapsulated, /inches,xsize=xs_ps, ysize=ys_ps
  endif
  if (j ge 9 AND j lt 99) then begin
    device,filename=strcompress(outpath+"rhox0" + string(j+1) + ".eps", /remove_all), $
           /color, /encapsulated, /inches,xsize=xs_ps, ysize=ys_ps
  endif
  if (j ge 99) then begin
    device,filename=strcompress(outpath+"rhox" + string(j+1) + ".eps", /remove_all), $
           /color, /encapsulated, /inches, xsize=xs_ps, ysize=ys_ps
  endif

  openr, 3, data_dir+"rho_xy.dat"
  tem=fltarr(NY)
  for i=0,num_p-1 do begin
    for k=0,NX-1 do begin
      readu,3,tem
      rho_xy(i,k,0:NY-1)=tem
    endfor
  endfor
  close,3

  rhox_vec=fltarr(NX)
  for l=0, NX-1 do begin
    rhox_vec[l]=0 
    for i=0, NY-1 do begin
      rhox_vec[l]+=rho_xy[j,l,i]
    endfor  
  endfor
  plot, x*100, 1.0e4*rhox_vec(*), /xstyle, /ystyle, yrange=[0.00, 0.7], $
        ;xrange=[-pipe_radius*100, pipe_radius*100], xtitle="x [cm]", $
        ytitle="density [arb. units]"
ENDFOR
ENDIF ; plot xprofile?


IF 1 THEN BEGIN ; plot xxs?
  device ,filename=strcompress(outpath+"xxs.eps", /remove_all), /color, /encapsulated, /inches, $
          xsize=xs_ps, ysize=ys_ps
  rho_xxs=fltarr(num_p, NX, NX)

  for j=0,num_p-1 do BEGIN
    openr, 3, data_dir+"xxs.dat"
    tem=fltarr(NX)
    for i=0,num_p-1 do begin
      for k=0,NX-1 do begin
        readu, 3, tem
        rho_xxs(i,k,0:NX-1)=tem
      endfor
    endfor
    close,3
  endfor
                        
  plot, pics(0,0,0:numprocs*Np-1)*100, pics(j,1,0:numprocs*Np-1)*1000.0, $
        /ystyle, psym=3, symsize=0.1, $
        xrange=[-pipe_radius*100,pipe_radius*100], yrange=[-xsmax,xsmax],$
   		;xrange=[-10,10], yrange=[-10,10],xstyle=1,$
    	xtitle="x [cm]", ytitle="x' [mrad]", /nodata
  		vline, 7
  arr = fltarr(NX,NX)
  arr(*,*) = rho_xxs(0,*,*)
  rho_xxs(0,*,*) = smooth(arr, smoothing, /edge_truncate)
  picmax = max(abs(rho_xxs(0,*,*)))
  picmin = picmax*1.0e-1      
  

  ang=2*!pi*indgen(360,/float)/359   
  gammax=(1+alpx^2)/betax  
  x0=x0-amp0  
  x0prime=x0prime-ampp0  
  e=gammax*x0^2+2*alpx*x0*x0prime+betax*x0prime^2

  
  t=findgen(cells) 
  ;inj_phase=acos(x0/sqrt(betax*e)) 
  inj_phase=-asin(sqrt(betax/e)*x0prime+alpx*x0/sqrt(betax*e))
  x2=sqrt(betax*e)*cos(inj_phase+2*!pi*Q_hor*t) 
  y2=sqrt(betax*e)*sin(inj_phase+2*!pi*Q_hor*t)
  z2=-(y2+alpx*x2)/betax   
 
  m=50
  time=findgen(m)/100  
  a0=findgen(m)
  ap0=findgen(m)
  a0[0]=amp0 
  ap0[0]=ampp0
  n=1
  for j=1, m-1 do begin
    if (j lt 5) then begin
      a0[j]=amp0-delAmp*n/4
      ap0[j]=ampp0-delAmp*ampp0/amp0*n/4
      n=n+1    
    endif
    if (j gt 4) then begin
      a0[j]=amp0-delAmp
	  ap0[j]=ampp0-delAmp*ampp0/amp0
    endif	
  endfor 
  x4=sqrt(betax*e)*cos(inj_phase+2*!pi*Q_hor*time) 
  y4=sqrt(betax*e)*sin(inj_phase+2*!pi*Q_hor*time)
  z4=-(y4+alpx*x4)/betax                              
   
  for j=0,num_p-1 do begin
  contour, rho_xxs(j,*,*), x*100, xs, /xstyle, /ystyle, $
            nlevels=clevels, /fill, $
           levels=picmin+findgen(clevels+1)*(picmax-picmin)/clevels, $
           C_COLORS=20+findgen(clevels+1)*(254-20)/clevels, /overplot
  endfor   
  x2[0]=x2[0]+amp0
  z2[0]=z2[0]+ampp0
  oplot, (x2)/1e-2,  (z2)/1e-3, psym=2   
  
  ;oplot, (x4+a0)/1e-2,  (z4+ap0)/1e-3, psym=0



  


ENDIF                           ; plot xxs?

device, /close
  spawn, 'display '+outpath+"xxs.eps"+' &'        
end
