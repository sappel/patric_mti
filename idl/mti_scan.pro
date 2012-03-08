pro mti_scan  
	
   	!P.CHARSIZE=1.7
	!P.CHARTHICK=1.0
	!P.THICK=4.0
	!Y.THICK=4.0
	!X.THICK=4.0    
	LoadCT, 13 
	
    tvlct, 255, 0, 0, 10 ;rot
	tvlct, 0, 0, 255, 11 ;blau
	tvlct, 0, 255, 0, 12 ;gruen
	
    !p.multi=0;[0,2,2]                
   
	
	;path='/d/bhs01/appel/patric/eickhoff/film_1_sp/47/'
	path='/d/bhs01/appel/patric/proton_sc/'      
	fileName='mti.dat'
	openr, lun, path+fileName, /get_lun
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
	turn=findgen(pmax)
	int=findgen(pmax)  
	loss=findgen(pmax)
	for i=0L, pmax-1, 1 do begin
	  turn[i] = data[3*i]
	  int[i] = data[1+3*i] 
	  loss[i] = data[2+3*i]
	endfor  
	
	;path='/d/bhs01/appel/patric/eickhoff/film_1_nsp/47/' 
	path='/d/bhs01/appel/patric/proton_sc2/'   
	fileName='mti.dat'
	openr, lun, path+fileName, /get_lun
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
  
    Ni = int[0]      
	entry_device=!d.name
	set_plot,"ps"
	!P.FONT=0   
	DEVICE,/TIMES     
    fileps='/u/sappel/tmp/scan.eps'                                 
	device,filename=fileps,/color,/encapsulated, xsize=20, ysize=15 

	    plot, turn, int, psym=0, xtitle='turns', ytitle=textoidl('N'), xmargin=[7,2.2], ymargin=[3.0,0.5], xrange=[0,100], /nodata
	    oplot, turn, int, psym=-1, color=10
		oplot, turn, int2, psym=-2, color=11      
		;oplot, turn, turn   
		;vline, 10
        ;plot, turn, loss, psym=-1, yrange=[0,1], $
		;xtitle='turns', ytitle=textoidl('loss'), xmargin=[7,2.2], ymargin=[3.0,0.5] 
		;plot, turn2, int2/1.e10, psym=-1, $
		;xtitle='turns', ytitle=textoidl('N/10^{10}'), xmargin=[7,2.2], ymargin=[3.0,0.5]
        ;plot, turn2, 1-loss2, psym=-1, yrange=[0,1], $
		;xtitle='turns', ytitle=textoidl('loss'), xmargin=[7,2.2], ymargin=[3.0,0.5]    
		!P.CHARSIZE=2.2
		legend, [ "sp", "sp2"], $
		psym=[-1,-2], color=[10,11], box=1, /bottom, /center

    device,/close
	set_plot, entry_device  
    
	!p.multi=0
	spawn, 'gv '+fileps+' &'
end	



