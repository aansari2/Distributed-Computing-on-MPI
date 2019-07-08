subroutine linspace(a,b,npts,x)
    integer, intent(in) :: npts
	real(8), intent(in) :: a, b
	real(8), intent(inout), dimension(npts) :: x
    integer :: i
    do i = 1,npts
        x(i) = a+(dble(i)-1.0D0)/(dble(npts)-1.0D0)*(b-a)
    end do
    x(npts) = b
end subroutine linspace

function func(x, y)
	real(8) :: x, y, func
    func = exp(x+y)*((x*x+3.0D0*x)*(y*y-y)+(y*y+3.0D0*y)*(x*x-x))
    return
end function

function sol(x, y)
    real(8) :: x, y, sol
    sol = exp(x+y)*(x*x-x)*(y*y-y)
    return
end function

program jacobicoarray
	real(8):: u(66,66)[*], f(66,66)[*], v(66,66)[*], s(66,66)[*], b(66,66)[*],&
	lerror[*], diff[*], gerror, gdiff
	real(8):: x(66)[*] !!
	real(8):: y(66)[*] !!
	real(8), CODIMENSION[*] :: xmin, xmax, ymin, ymax, threshold
	integer, CODIMENSION[*] :: rank, p, q, i, j, xlim, ylim, ix, iterCount
	character(len=32) :: buffer[*]
	buffer = ""

	call get_command_argument(1, buffer)
	read(buffer, *) iterCount
	call get_command_argument(2, buffer)
	read(buffer, *) threshold

	rank = this_image()

	! Grid Location
	p = mod(rank-1,4)
	q = (rank-1)/4

	! Grid Limits
	xmin = dble(p)/4.0D0
	xmax = 0.25D0 + 1.0D0/256.0D0 + dble(p)/4.0D0 !!
	ymin = dble(q)/4.0
	ymax = 0.25D0 + 1.0D0/256.0D0 + dble(q)/4.0D0 !!

	! Define Local Grid
	call linspace(xmin, xmax, 66, x)!!
	call linspace(ymin, ymax, 66, y)!!



	! Define f
	do i = 1,66!!
		do j = 1,66!!
			u(i,j)[rank] = 0.0*dble(rank)!!
			v(i,j)[rank] = 0.0*dble(rank)!!
			f(i,j)[rank] = func(x(i),y(j))
			s(i,j)[rank] = sol (x(i),y(j))
		end do
	end do
	u=0.0
	v=0.0
	b=0.0
	k = 0
	do while (1)
		if(mod(k,10)==0) then
			if (k/=0) then
				! ######      ###    ######## ##     ## ######## ########        ######## ########  ########      
				!##    ##    ## ##      ##    ##     ## ##       ##     ##       ##       ##     ## ##     ##     
				!##         ##   ##     ##    ##     ## ##       ##     ##       ##       ##     ## ##     ##     
				!##   #### ##     ##    ##    ######### ######   ########        ######   ########  ########      
				!##    ##  #########    ##    ##     ## ##       ##   ##         ##       ##   ##   ##   ##       
				!##    ##  ##     ##    ##    ##     ## ##       ##    ##        ##       ##    ##  ##    ##  ## 
				! ######   ##     ##    ##    ##     ## ######## ##     ##       ######## ##     ## ##     ## ors 
				
				!find diff Exact minus numerical
				diff = maxval(abs(s(2:65,2:65)-u(2:65,2:65)))
				! All reduce function from 'http://www.training.prace-ri.eu/uploads/tx_pracetmo/L04_Experiences.pdf'
				! as comax function wasn't working
				!all reduce to gdiff
				sync all
				gdiff = diff
				do i=1,num_images()
					gdiff=max(gdiff,diff[i])
				end do
				!find error b - u
				lerror = maxval(abs(b(2:65,2:65)-u(2:65,2:65)))
				!All reduce to gerror
				sync all
				gerror = lerror
				do i=1,num_images()
					gerror=max(gerror,lerror[i])
				end do
				!set b = u
				b = u
				
				!if (rank==1) then
				!	write(*,*) k,gerror,gdiff
				!end if
				sync all !Synchronize
				!if statement --> print and break
				if ((gerror<threshold).AND.(k>=iterCount)) then
					if (rank==1) then
						write(*,*) 'd = ',gerror,'	c = ',gdiff
					end if
					exit
				end if
			end if
		end if
		k = k + 1
		!########     ###    ########    ###          ######## ##     ##  ######  ##     ## 
		!##     ##   ## ##      ##      ## ##         ##        ##   ##  ##    ## ##     ##       
		!##     ##  ##   ##     ##     ##   ##        ##         ## ##   ##       ##     ##       
		!##     ## ##     ##    ##    ##     ##       ######      ###    ##       #########   
		!##     ## #########    ##    #########       ##         ## ##   ##       ##     ##       
		!##     ## ##     ##    ##    ##     ##       ##        ##   ##  ##    ## ##     ## ###       
		!########  ##     ##    ##    ##     ##       ######## ##     ##  ######  ##     ##ange
		sync all
		if (p /= 3) then
			u(66,1:66)[rank] = u(2,1:66)[rank+1]
		end if

		if (p /= 0) then
			u(1,1:66)[rank] = u(65,1:66)[rank-1]
		end if	

		if (q /= 3) then
			u(1:66,66)[rank] = u(1:66,2)[rank+4]
		end if

		if (q /= 0) then
			u(1:66,1)[rank] = u(1:66,65)[rank-4]
		end if
		sync all
		!      ##    ###     ######   #######  ########  ####       #### ######## ######## ########  
		!      ##   ## ##   ##    ## ##     ## ##     ##  ##         ##     ##    ##       ##     ## 
		!      ##  ##   ##  ##       ##     ## ##     ##  ##         ##     ##    ##       ##     ## 
		!      ## ##     ## ##       ##     ## ########   ##         ##     ##    ######   ########  
		!##    ## ######### ##       ##     ## ##     ##  ##         ##     ##    ##       ##   ##   
		!##    ## ##     ## ##    ## ##     ## ##     ##  ##         ##     ##    ##       ##    ##   ###
		! ######  ##     ##  ######   #######  ########  ####       ####    ##    ######## ##     ## ation
		xlim = merge(64,65,p==3) !!
		ylim = merge(64,65,q==3) !!
		
		v(2:xlim,2:ylim) = 0.25D0*(u(3:xlim+1,2:ylim) + u(1:xlim-1,2:ylim) + u(2:xlim,3:ylim+1) + &
		u(2:xlim,1:ylim-1) - f(2:xlim,2:ylim)/65536.0D0) !!
		
		
		u = v
	end do
end program jacobicoarray

! Author Adil Ansari
!ssso++++++ooooyo+ooo++++/++////oososo/:::::::::--::--------------------::::::-::--:::::
!ooooooo++ooo++++++o+++///+++++/+ossyyso+//:://::::::::::----------------:::::::::::::::
!oooooo+++oo++++//:::/:/++////+/++++ydhy+ooosssssyyssyhdy+/:-------------:::::::::::::::
!ooo++o+++oo++++++//////+o+//++++///ohddssyhdmmNNNmdddddmmdysss+/::------:::::::::::::::
!oo++//++oo+oooo++++/////+//++/+o+//ydNNNmmNNNmmmmNmmmmmmNNNNNNmh+:::::::::::::::::::://
!+++///////////++osso++/////++++oooyNMMMMNNNNNmmmNNNmmmmNMMMMMMNNmo/::::::::::::::::////
!+ooo+++++o////+ossyyyo+/+/////+osymMMMMMMNNNNNMMMMNNNNNNNNMMMMMMMmyo:::::::::::::::////
!sssoooo++o///++oso+shysoooo+++osydNMMMMMMMMMMMMMMNmmmmmNmmmNMMMMMMNmh+::::::::::///////
!soooooo+////++/:/+++++osysyyhhhhdmNMMMMNMMMMMMMMMNNNNNNNNNNmMMMMMMMMNd/::/:::://///////
!++++++++///////++++//+++++oydNNmmmNNmdddmmmmmmNmmdNNMMMMMMNNMMMMMMMMMNs////////////////
!+++++/++++++/://++/////////ymMNNmddhs+//oooooo++++ydmNMMMMNNNNMMMMMMMMNs+////////////++
!+++/::///////////////////+ohNMNmmmmh+:-.---..``.-:+yhdmNNNNNNNMMMMMMMMMho//////////++++
!o+///////////:/++++++//::-+hNNNNNNNy/-.````    ```.:+oyhdmmNNNNMMMMMMMds/:://////++++++
!+///:---:/:::/++////o+++/:odNNMNNmd/:..`````````````..-:+oyhdmmNNMMMMMo/::::///++++++++
!++++/::::/:-:://///+oosso/sdNNMNho/----....`````````````.-:+yhddmMMMMN+//+/+/:://++oo++
!//////+////////:::/+//::::ohNMNy+:-://+osssoo/:--...-...---:+syhdMMMMm///+ss+/---:::///
!//////ooooo++/:-:::://:::/+ymMhyyso+++oydmmmmmdyo--+ossyyyssosyyhNMMNy:://o++/::::--:::
!///+/+sssoo+/:::::/+++///ssyhh+hdh/oyhhdmmmmmmmmmyddddmNNNmmdmmdhNNNdo::://++/::/::::--
!/++++ooooooo++///++++++///:+s+-+ss--/+ossyyyyyhyo/hdddmNNNNmmdddNmmhso/////+////////:--
!+ooo++++///+oo/+//+++++//..:+:.-+o:...-:::---+/-``/ssossyyhdmdhhmdmdy+++++++/////:://--
!ooooooo+/////++//+osoo++/...::...-:......--::.`` `./+/--::/+oooo+hmdysssoos+/://::://::
!//++ooo+//://///://++++//.``.:-........-/+:-```  `.-/+/:--:://::/hdhsosooos+//+o////+//
!//++//+++++++//o+////+///-.`.::--.....-/o:--.....-::/+o/::::://+odyo//:::////++o//+/://
!///+//oo+/////://////++++/:-:///::::-:/+--/sysoosyhhhs+::://+oosss+/::::::://+++::///++
!:///+so/:://///+////////://////////:://:-.-:/yhddmmddhs///++sssss+///:-:/:-:::::-:/+s+/
!///++///+//////+++++//:::::::++/////:/+syo+///osyyhhddy/+osyyys+::::://+//::::::::/+y+/
!/+++///+++////+ooooo++////:::+++///////://::::/+oydmNmdoosyyyyo/:/++++++++/:::---:++o/:
!oo++/////:///++soooso+++/////+oo+++//::--------:/ohhdhysyyhhyo+//:/ossshdhy//::----::-:
!ooo++/::::////++///++/////////ssso+/::---::+oyyyysoooooshhhhs++/:-:/osshdho/:---:---:-:
!////////:://////////////:::::-oyyyo+/:--..---::::/::/+ohhhs+ohhy/:/++/+++/::::::::-----
!::::-:::::-:///+++/+ooo++/:--.-oyhhys+/--.....---::+sydds+//+oss+/////:///:::://::::::-
!::://:::::::::/+++++:--..---...-+ydddhyo++//+++oooyhdddy+/:////////++//+++/:///+o++////
!`.-://///++++++ooo+-`````.......-+yddddddddddmmmmmmmmddhhyo:::://:::::/+o+/://///:::::/
!:::/+oooosooooosss+.``````.......-:shhdddmmmmmmmmmdddddmNNmy//////////:::::////::::::::
!sssssssssoooooossso:`````.........-/oyhhdddmmmmmdddddddmNNNNmdyo///////::::+++////:::::
!ssssssssssoooossssso.```````......--/oydddmmmdddddddddmNNMMMMMNNdo+////:::::////::--:::
!ssssssssssossssoossso:.``````.....-:/+shddddddddddddddNMMMMMMNNmdhso++/:-------::-----:
!ssssssssssosssssoosssso/-````......--:/+hdddddddddmmNMMMMMMNmhysooooo++/:---.......-..-
!sssssssssssssssssssssssss+/:--...----.../shdmmmmNNMMMMNMNNmhsoooooooooo++:---....----..
!ssssssssssssssssssssssssyyyysso++//:::::/syhmNNMMMMMMNNmdysooooooooooo+++++:----:///:--