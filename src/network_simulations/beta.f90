subroutine network_sim(ttot,nlaps,length_ts,Nff_act,spatial_neurons,spatial_neurons_ec,&
    &Ncontext_act,spatial_matrix,nonspatial_matrix,aei,aii,aie,matr_cont_inh,ceffinh,base_dir,Nexc,Nff,Ncontext,Ninh,&
    &output_matrix,output_matrix_inh)
    implicit none
    INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
    real(kind=DP), intent(in) :: ttot
    CHARACTER(len=64), intent(in) :: base_dir
    integer,intent(in) :: Nexc,Ninh,Ncontext,Nff,nlaps,length_ts,Nff_act,spatial_neurons,Ncontext_act,spatial_neurons_ec
    integer(kind=4), dimension(Nexc*8,Nff), intent(in):: spatial_matrix
    integer, dimension(Nexc*8,Ncontext), intent(in):: nonspatial_matrix
    integer, dimension(Ninh,Ninh), intent(in) :: aii
    integer, dimension(Nexc,Ninh), intent(in) :: aei
    integer, dimension(Ninh,Nexc), intent(in) :: aie
    integer, dimension(Ninh,Nff), intent(in) :: ceffinh
    integer, dimension(Ninh,Ncontext), intent(in) :: matr_cont_inh
    real(kind=DP)::dt,tau,nran,geff,geffk,theta,iext,taus,v_reset,rb,rphi,phi,phii,tauinh,tausinh,ginhffk
    real(kind=DP),dimension(Nexc) :: vexc,Isynff,Isynei,Iavgff,Iavgcont,Iavginh,synexc,synexccont,syneinh_old,syneinh
    real(kind=DP),dimension(Ninh) :: vinh,Iavgcontinh,Isyni,Isynie,synexcinh,synexccontinh,syninhe,syniinh,syniinh_old,&
    &syninhe_old
    real(kind=DP),dimension(Nff) :: ratesff
    real(kind=DP),dimension(Nff) :: ratescontext
    integer(kind=4),dimension(Nexc,Nff) :: ceff
    integer(kind=4),dimension(Nexc,Ninh) :: conie
    integer(kind=4),dimension(Nexc,Ncontext) :: matr_cont
    integer(kind=4),dimension(Nff):: neur_ind
    integer(kind=4),dimension(Ncontext)::neur_ind_cont
    integer(kind=4),dimension(Ninh,Ninh):: conii
    integer(kind=4),dimension(Ninh,Nexc)::conei
    integer(kind=4),dimension(Nff)::ncount,ncountfinh
    integer(kind=4),dimension(Ncontext)::ncountcont,ncountcontinh
    integer(kind=4),dimension(Nff,Nexc)::conff
    integer(kind=4),dimension(Nff,Ninh)::conffinh
    integer(kind=4),dimension(Ncontext,Nexc)::ccont
    integer(kind=4),dimension(Ncontext,Ninh)::ccontinh
    integer(kind=4),dimension(Nexc)::ncountexc
    integer(kind=4),dimension(Ninh)::ncountinh,ncountinhinh
    integer(kind=4),dimension(Nff_act)::active_ca3
    integer(kind=4),dimension(spatial_neurons):: tuned_neurs
    integer(kind=4),dimension(spatial_neurons_ec):: tuned_neurs_ec
    real(kind=DP),dimension(spatial_neurons)::phivector
    real(kind=DP),dimension(spatial_neurons_ec)::phivector_ec
    logical,dimension(Nff)::is_spatial
    logical,dimension(Ncontext)::is_spatial_ec
    integer(kind=4),dimension(Ncontext_act)::active_context    
    integer:: Nstep,i,j,l,nl,avg_steps,Kef,Kii,Kie,Kei,nz,nr,Kefinh,tt
    real(kind=DP)::pi=4.d0*datan(1.d0),period,freq,muc,sigmac   
    real(kind=DP)::geik,gei,giik,gii,giek,gie,t_avg,start,finish,starttot,dtp,Iextinh,ginhff,sigma_bin,max_rate,kvm
    real(kind=DP)::gcont,gcontinh,meanord,sigmaord,gcontk,gcontinhk,pos_x,length_track,velocity_x,x_i
    integer::length_inh,length_exc,looper,kj,n_t_kj,nsess=0,it_spat
    CHARACTER(len=32) :: file_spat,file_nonspat
    integer::neur_loc_cont,spike_bin,spike_exc
    real,dimension(length_ts,Nexc),intent(out) :: output_matrix
    real,dimension(2*length_ts,Ninh),intent(out) :: output_matrix_inh
    integer,dimension(Nexc) :: dk
    integer,dimension(Ninh) :: dinh
    integer,dimension(20) :: selected_cells,selected_cells_inh
    
    !call init_random_seed()
    
    !Nff_act=nint(ca3_spars*Nff)
    !Ncontext_act=nint(cont_spars*Ncontext)
    !spatial_neurons=nint(frac_tun*Nff_act)
    !ttot=1120.d3 !!in milliseconds
    !ttot=5.d3
    dt=0.1d0
    Nstep=nint(ttot/dt)
    !Nstep=10
    Kef=500 !!average number of connections per neuron
    Kii=500
    Kie=500
    Kei=500
    Kefinh=500
    meanord=500.d0
    sigmaord=47.d0
    geffk=3.5d0
    geik=2.8d0
    giek=2.d0
    giik=2.d0
    ginhffk=3.d0
    gcontk=2.52744d0
    gcontinhk=2.d0
    !!carefoul
    geff=geffk/(dble(Kef)**0.5)
    gei=geik/(dble(Kei)**0.5)
    gie=giek/(dble(Kie)**0.5)
    gii=giik/(dble(Kii)**0.5)
    gcont=gcontk/(dble(meanord)**0.5)
    gcontinh=gcontinhk/(dble(meanord)**0.5)
    ginhff=ginhffk/(dble(Kefinh)**0.5)
    theta=10.d0
    iext=0.d0
    tau=10.d0
    tauinh=10.d0
    taus=5.d0
    tausinh=5.d0
    v_reset=-0.d0
    rb=5.d0
    rphi=40.d0
    t_avg=20.d0 !milliseconds
    avg_steps=nint(t_avg/dt)
    dtp=0.d0
    muc=3.d0
    sigmac=0.77d0
    dk=1
    dinh=1
    conff=0
    conffinh=0
    ncount=1
    ncountcont=1
    ncountcontinh=1
    ncountfinh=1
    synexc=0.d0
    synexcinh=0.d0
    synexccont=0.d0
    synexccontinh=0.d0
    conie=0
    ncountexc=1
    syninhe=0.d0
    conei=0
    ncountinh=1
    syneinh=0.d0
    syneinh_old=0.d0
    conii=0
    ncountinhinh=1
    syniinh=0.d0
    syniinh_old=0.d0    
    Iavgff=0.d0
    Iavginh=0.d0
    Isynff=0.d0
    Isynei=0.d0
    Isyni=0.d0
    Isynie=0.d0
    vexc=0.1d0
    vinh=0.d0
    period=7.d3
    freq=1.d0/period
    Iavgcont=0.d0
    Iavgcontinh=0.d0
    ratescontext=13.27d0
    matr_cont=0
    velocity_x=12.d-3
    length_track=84.d0 !! in centimeters
    sigma_bin=3.d0
    max_rate=90.d0
    kvm=length_track**2/(sigma_bin**2*4.d0*pi**2)
    !nlaps=20 
        !!!non spatial variance calculation
    sigmaord=20.916500663d0
    !!initial conditions!!!
    output_matrix=0.
    output_matrix_inh=0.
    sigmaord=20.916500663d0
    !!initial conditions!!!
    do i=1,Nexc
        call random_number(nran)
        vexc(i)=nran*5.d0
    end do
    
    do i=1,Ninh
        call random_number(nran)
        vinh(i)=nran*5.d0
    end do
    !!list of neuron indeces!!!
    
    do i=1,Nff
        neur_ind(i)=i
    end do
    
    do i=1,Ncontext
        neur_ind_cont(i)=i
    end do

    !!cells with saved currents, voltages, etc
    do i=1,20
        call random_number(nran)
        selected_cells(i)=nint(nran*Nexc)
    end do
    
    !!cells with saved currents, voltages, etc (inhibitory)
    do i=1,20
        call random_number(nran)
        selected_cells_inh(i)=nint(nran*Ninh)
    end do 
 
    open(unit=89, file=trim(base_dir) // '/data/outputs_simulation/indices_selected.dat',action='write')
    write(89,*) (selected_cells(j),j=1,20)
    close(89)
    
    open(unit=92, file=trim(base_dir) // '/data/outputs_simulation/indices_selected_inh.dat',action='write')
    write(92,*) (selected_cells_inh(j),j=1,20)
    close(92)

!    !!!!ca3-I connectivity matrix 
!    open(unit=188,file="ceffinh.dat",action='read')
!    do tt=1,Ninh
!        read(188,*) ceffinh(tt,:)
!    end do 
!    close(188)
    
    
    do i=1,Ninh
        do j=1,Nff
            !call random_number(nran)
            if (ceffinh(i,j)==1) then 
                !ceffinh(i,j)=1
                conffinh(j,ncountfinh(j))=i
                ncountfinh(j)=ncountfinh(j)+1
            end if
        end do 
    end do
        
    
    
!   !!!!contextual to inhibitory connectivity matrix 
!    open(unit=188,file="matr_cont_inh.dat",action='read')
!    do tt=1,Ninh
!        read(188,*) matr_cont_inh(tt,:)
!    end do 
!    close(188)
!!    
    
    do i=1,Ninh
        do j=1,Ncontext
            !call random_number(nran)
            if (matr_cont_inh(i,j)==1) then 
                !ceffinh(i,j)=1
                ccontinh(j,ncountcontinh(j))=i
                ncountcontinh(j)=ncountcontinh(j)+1
            end if
        end do 
    end do
        
    
    

    
!    
        !!!connectivity matrix aie!!!
    do i=1,Ninh
        do j=1,Nexc
            !call random_number(nran)
            if (aie(i,j)==1) then 
                conie(j,ncountexc(j))=i
                ncountexc(j)=ncountexc(j)+1
            end if
        end do 
    end do
    
        !!!connectivity matrix aii!!!
    !aii=0
    do i=1,Ninh
        do j=1,Ninh
            !call random_number(nran)
            if (aii(i,j)==1) then 
                !aii(i,j)=1
                !aii(j,i)=1
                conii(j,ncountinhinh(j))=i
                ncountinhinh(j)=ncountinhinh(j)+1
            !else
                !aii(i,j)=0
                !aii(j,i)=0
            end if
        end do 
    end do

        !!!connectivity matrix aei!!!
    do i=1,Nexc
        do j=1,Ninh
            !call random_number(nran)
            if (aei(i,j)==1) then 
                conei(j,ncountinh(j))=i
                ncountinh(j)=ncountinh(j)+1
            
            end if
        end do 
    end do
    ratesff=10.d0 !poisson rates
    !!assign a phase to each input
    phi=0.d0
    do i=1,Nff
        phii=-pi+(i-1)*pi*2.d0/dble(Nff)
        ratesff(i)=rb+rphi*(1.d0+dcos(phi-phii))
        ratesff(i)=40.d0
    end do 

    open(unit=60,file=trim(base_dir) // '/data/outputs_simulation/membrane_potentials.dat',action='write')
    open(unit=70,file=trim(base_dir) // '/data/outputs_simulation/input_currents_cont_e.dat',action='write')
    open(unit=71,file=trim(base_dir) // '/data/outputs_simulation/input_currents_spat_e.dat',action='write')
    open(unit=72,file=trim(base_dir) // '/data/outputs_simulation/input_currents_inh_e.dat',action='write')
    open(unit=73,file=trim(base_dir) // '/data/outputs_simulation/membrane_potentials_inh.dat',action='write')
    open(unit=74,file=trim(base_dir) // '/data/outputs_simulation/input_currents_cont_i.dat',action='write')
    open(unit=75,file=trim(base_dir) // '/data/outputs_simulation/input_currents_spat_i.dat',action='write')
    open(unit=76,file=trim(base_dir) // '/data/outputs_simulation/input_currents_exc_i.dat',action='write')
    open(unit=77,file=trim(base_dir) // '/data/outputs_simulation/input_currents_inh_i.dat',action='write')
    
    Iext=0.d0
    Iextinh=0.d0
    call cpu_time(start)
    !!!!assing contextual inputs
    ratescontext(:)=0.d0
    call randsampl(neur_ind_cont,active_context,Ncontext,Ncontext_act)
    ratescontext(active_context)=13.27d0
    open(unit=1200,file=trim(base_dir) // '/data/outputs_simulation/rates_EC_one_dir.dat',action='write')
    do nl=1,Ncontext_act
        write(1200,*) active_context(nl)
    end do 
    close(1200)
    
    !!!!select subset of ECIII cells with spatial tuning!!!!
    call randsampl(active_context,tuned_neurs_ec,Ncontext_act,spatial_neurons_ec)
    open(unit=1121,file=trim(base_dir) // '/data/outputs_simulation/tuned_ec3_cells.dat',action='write')
    do nl=1,spatial_neurons_ec
        write(1121,*) tuned_neurs_ec(nl)
    end do 
    close(1121)
    open(unit=8661,file=trim(base_dir) // '/data/outputs_simulation/input_phase_ec.dat',action='write')
    do nl=1,spatial_neurons_ec
        !call random_number(nran)
        phivector_ec(nl)=-pi+(nl-1)*pi*2.d0/dble(spatial_neurons_ec)
        write(8661,*) phivector_ec(nl)
    end do
    close(8661)
    
    
    is_spatial_ec(:)=.False.
    do nl=1,Ncontext
        if(any(tuned_neurs_ec==nl)) then
                is_spatial_ec(nl)=.True.
        end if
    end do  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!assign the input firing rates
    ratesff(:)=0.d0
    call randsampl(neur_ind,active_ca3,Nff,Nff_act) !!neurons active but with no tuning
    open(unit=1100,file=trim(base_dir) // '/data/outputs_simulation/rates_ca3_one_dir.dat',action='write')
    do nl=1,Nff_act
        write(1100,*) active_ca3(nl)
    end do 
    close(1100)
    call randsampl(active_ca3,tuned_neurs,Nff_act,spatial_neurons) !!neurons active with  tuning
    open(unit=1110,file=trim(base_dir) // '/data/outputs_simulation/tuned_ca3_cells.dat',action='write')
    do nl=1,spatial_neurons
        write(1110,*) tuned_neurs(nl)
    end do 
    close(1110)
    
    
    ratesff(active_ca3)=13.27d0
    open(unit=868,file=trim(base_dir) // '/data/outputs_simulation/input_phase.dat',action='write')
    do nl=1,spatial_neurons
        !call random_number(nran)
        phivector(nl)=-pi+(nl-1)*pi*2.d0/dble(spatial_neurons)
        write(868,*) phivector(nl)
    end do
    close(868)
    
    is_spatial(:)=.False.
    do nl=1,Nff
        if(any(tuned_neurs==nl)) then
                is_spatial(nl)=.True.
        end if
    end do  
    
    do i=1,Nstep
        call cpu_time(starttot)
        pos_x=-length_track/2.d0+mod(velocity_x*i*dt,length_track)
        phi=2.d0*pi*pos_x/(length_track)
        
        !spikes_ff=0
        synexc=synexc-synexc*dt/taus
        synexcinh=synexcinh-synexcinh*dt/taus
        synexccont=synexccont-synexccont*dt/taus
        synexccontinh=synexccontinh-synexccontinh*dt/taus
        it_spat=1
        do nl=1,Nff
            if (is_spatial(nl)) then
                phii=phivector(it_spat)
                x_i=phii*length_track/(2.d0*pi)
                ratesff(nl)=max_rate*dexp(kvm*dcos(phi-phii))/(dexp(kvm))+rb
                it_spat=it_spat+1
            end if 

            call random_number(nran)
            if (ratesff(nl)*dt/(1000.d0)>=nran) then 
                do nz=1,ncount(nl)-1
                    if (conff(nl,nz)==0) then 
                        print*,"error: unassigend postsynaptic cell",nl,nz
                    end if 
                    synexc(conff(nl,nz))=synexc(conff(nl,nz))+1.d0/taus
                    
                end do
                do nr=1,ncountfinh(nl)-1
                    if (conffinh(nl,nr)==0) then 
                        print*,"error: unassigend postsynaptic cell",nr,nz
                    end if 
                    synexcinh(conffinh(nl,nr))=synexcinh(conffinh(nl,nr))+1.d0/taus
                end do
            end if
        end do
        
        it_spat=1
        do nl=1,Ncontext_act
            neur_loc_cont=active_context(nl)

            if (is_spatial_ec(neur_loc_cont)) then 
                phii=phivector_ec(it_spat)
                x_i=phii*length_track/(2.d0*pi)
                ratescontext(neur_loc_cont)=max_rate*dexp(kvm*dcos(phi-phii))/(dexp(kvm))+rb
                it_spat=it_spat+1
            end if 
            call random_number(nran)
            if (ratescontext(neur_loc_cont)*dt/(1000.d0)>=nran) then 
                do nz=1,ncountcont(neur_loc_cont)-1
                    if (ccont(neur_loc_cont,nz)==0) then 
                        print*,"error: unassigend postsynaptic cell",nl,nz
                    end if 
                    synexccont(ccont(neur_loc_cont,nz))=synexccont(ccont(neur_loc_cont,nz))+1.d0/taus
                end do
                do nr=1,ncountcontinh(neur_loc_cont)-1
                    if (ccontinh(neur_loc_cont,nr)==0) then 
                        print*,"error: unassigend postsynaptic cell",nr,nz
                    end if 
                    synexccontinh(ccontinh(neur_loc_cont,nr))=synexccontinh(ccontinh(neur_loc_cont,nr))+1.d0/taus
                end do
               
            end if
        end do
     
!!the degree of each neuron is a wigthed sum of the previous degree and the new one!!!! 
        if (mod(i-1,nlaps*nint(period/dt))==0 .and. i*dt<1120.d3) then
            nsess=nsess+1
            ncountcont=1

            ceff=spatial_matrix(Nexc*(nsess-1)+1:Nexc*nsess,:)
            matr_cont=nonspatial_matrix(Nexc*(nsess-1)+1:Nexc*nsess,:)
            
            

            
            do tt=1,Ncontext
                length_exc=sum(matr_cont(:,tt))
                length_inh=sum(matr_cont_inh(:,tt))
                ncountcont(tt)=length_exc+1
                ncountcontinh(tt)=length_inh+1
                looper=1
                do kj=1,Nexc
                    if (matr_cont(kj,tt)==1) then 
                        ccont(tt,looper)=kj
                        looper=looper+1
                    end if
                end do
                looper=1
                do kj=1,Ninh
                    if (matr_cont_inh(kj,tt)==1) then 
                        ccontinh(tt,looper)=kj
                        looper=looper+1
                    end if
                end do
              
            end do 
            
            
            !!!!!postsyn spatial!!!!
            n_t_kj=0
            do tt=1,Nff
                length_exc=sum(ceff(:,tt))
                ncount(tt)=length_exc+1
                looper=1
                do kj=1,Nexc
                    if (ceff(kj,tt)==1) then 
                        conff(tt,looper)=kj
                        looper=looper+1
                        if (kj==4000) then
                            n_t_kj=n_t_kj+1
                        end if 
                    end if
                end do
              
            end do 
            
            

    
        
        end if
        
        
        syninhe_old=syninhe
        syninhe=syninhe-syninhe*dt/taus
        spike_exc=0
        do l=1,Nexc
  
                Iavgff(l)=(Iavgff(l)+synexc(l))

                Iavginh(l)=(Iavginh(l)+syneinh_old(l))

                Iavgcont(l)=(Iavgcont(l)+synexccont(l))

            vexc(l)=lif(vexc(l),geff*synexc(l)-gei*syneinh(l)+gcont*synexccont(l),Iext,tau,dt)
            if (vexc(l)>theta) then
                spike_exc=spike_exc+1
                !write(50,*) l,pos_x,i*dt,phi
                output_matrix(dk(l),l)=i*dt
                !output_matrix(dk(l),l)=i
                vexc(l)=v_reset
                dk(l)=dk(l)+1
                do nz=1,ncountexc(l)-1
                    if (conie(l,nz)==0) then 
                        print*,"error: unassigend postsynaptic cell",nl,nz
                    end if 
                    syninhe(conie(l,nz))=syninhe(conie(l,nz))+1.d0/taus
                end do 
            end if 
        end do
        !!firing rate!!!!
        

        syneinh_old=syneinh
        syneinh=syneinh-syneinh*dt/tausinh
        !spike_bin_vec_exc(mod(i,avg_steps)+1)=dble(spike_exc)/4000/dt*1000
        syniinh_old=syniinh
        syniinh=syniinh-syniinh*dt/tausinh
        spike_bin=0
        do l=1,Ninh
            
            vinh(l)=lif(vinh(l),-gii*syniinh_old(l)+gie*syninhe_old(l)+ginhff*synexcinh(l)+gcontinh*synexccontinh(l)&
            &,Iextinh,tauinh,dt)
            Iavgcontinh(l)=(Iavgcontinh(l)+synexccontinh(l)/Nstep)
            if (vinh(l)>theta) then
                spike_bin=spike_bin+1
                !write(55,*) l,pos_x,i*dt,phi
                output_matrix_inh(dinh(l),l)=i*dt
                dinh(l)=dinh(l)+1
                vinh(l)=v_reset
               do nz=1,ncountinh(l)-1
                    if (conei(l,nz)==0) then 
                        print*,"error: unassigend postsynaptic cell",nl,nz
                    end if 
                    syneinh(conei(l,nz))=syneinh(conei(l,nz))+1.d0/tausinh
                end do
                do nz=1,ncountinhinh(l)-1
                    if (conii(l,nz)==0) then 
                        print*,"error: unassigend postsynaptic cell",nl,nz
                    end if 
                    syniinh(conii(l,nz))=syniinh(conii(l,nz))+1.d0/tausinh
                end do
            end if 
        end do

        
        !!mean field inhibitory synaptic current!!
!        if (mod(i,25)==0) then 
!            write(10,*) i*dt,sum(synexc)/Nexc*tau*geff,&
!            &sum(syneinh_old)/Nexc*tau*gei,sum(syninhe_old)/Ninh*tauinh*gie,&
!            &sum(synexcinh)/Ninh*tau*ginhff,sum(syniinh_old)/Ninh*tauinh*gii,&
!            &sum(synexccont)/Nexc*tau*gcont,sum(synexccontinh)/Ninh*tau*gcontinh,&
!            &synexc(1)*tau*geff,vexc(1),vinh(1),syneinh_old(1)*tau*gei,synexc(4000)



!        end if
        if (mod(i,100)==0) then
            write(60,*) (vexc(selected_cells(j)),j=1,20)
            write(70,*) (gcont*synexccont(selected_cells(j)),j=1,20)
            write(71,*) (geff*synexc(selected_cells(j)),j=1,20)
            write(72,*) (gei*syneinh(selected_cells(j)),j=1,20)
            write(73,*) (vinh(selected_cells_inh(j)),j=1,20)
            write(74,*) (gcontinh*synexccontinh(selected_cells_inh(j)),j=1,20)
            write(75,*) (ginhff*synexcinh(selected_cells_inh(j)),j=1,20)
            write(76,*) (gie*syninhe(selected_cells_inh(j)),j=1,20)
            write(77,*) (gii*syniinh(selected_cells_inh(j)),j=1,20)
        end if 

    end do
    call cpu_time(finish)

    close(60)
    close(70)
    close(10)
    close(20)
    close(50)
    close(3)
    close(4)
    close(22)
    close(24)
    close(73)
    close(55)
    close(71)
    close(72)
    close(73)
    close(74)
    close(75)
    close(76)
    close(77)
    return 
    contains
    
    real (kind=DP) function lif(v,isyn,iext,tau,dt) result (res)
    INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
    real(kind=DP)::v,isyn,iext,tau,dt
    res=v+(-v+Isyn*tau+iext)*dt/tau
    end function lif

    
    function normal(mean, sigma)

    ! mean  : mean of distribution
    ! sigma : number of standard deviations

    implicit none

    integer, parameter:: b8 = selected_real_kind(14)
    real(b8), parameter :: pi = 3.141592653589793239_b8
    real(b8) normal
    real rand_num
    real(b8) tmp
    real(b8) mean
    real(b8) sigma
    real(b8) fac
    real(b8) gsave
    real(b8) rsq
    real(b8) r1
    real(b8) r2
    integer flag
    save flag
    save gsave
    data flag /0/

    if (flag.eq.0) then
    rsq=2.0_b8
    do while(rsq.ge.1.0_b8.or.rsq.eq.0.0_b8) ! new from for do
      call random_number(rand_num)
      r1=2.0_b8*rand_num-1.0_b8
      call random_number(rand_num)
      r2=2.0_b8*rand_num-1.0_b8
      rsq=r1*r1+r2*r2
    enddo

    fac=sqrt(-2.0_b8*log(rsq)/rsq)
    gsave=r1*fac
    tmp=r2*fac
    flag=1
    else
    tmp=gsave
    flag=0
    endif

    normal=tmp*sigma+mean

    return
    endfunction normal
    
    real(kind=DP) function gaussian(x,sigma)
    real(kind=DP)::x,sigma
    real(kind=DP)::pi=4.d0*datan(1.d0)
    gaussian=(1.d0/(2.d0*pi*sigma**2)**0.5)*dexp(-(x)**2/(2.d0*sigma**2))
    return 
    endfunction gaussian 
    
    subroutine randsampl(x,a,n,k)
    integer,dimension(n)::x
    integer,dimension(k)::a
    integer::N,k,j,l,m
    real(kind=DP)::nran 
    m=0
    do j=1,n
        call random_number(nran)
        l=int((dble(n-j+1))*nran)+1
        if (l.gt. (k-m)) cycle
        m=m+1
        a(m)=x(j)
        if (m .ge. k) then

            exit
        end if 
    end do
    return
    end subroutine 
    
    real(kind=DP) function bessel_i0 (arg)

    !*****************************************************************************80
    !
    !! BESSEL_I0 evaluates the modified Bessel function I0(X).
    !
    !  Discussion:
    !
    !    The main computation evaluates slightly modified forms of
    !    minimax approximations generated by Blair and Edwards, Chalk
    !    River (Atomic Energy of Canada Limited) Report AECL-4928,
    !    October, 1974.  This transportable program is patterned after
    !    the machine dependent FUNPACK packet NATSI0, but cannot match
    !    that version for efficiency or accuracy.  This version uses
    !    rational functions that theoretically approximate I-SUB-0(X)
    !    to at least 18 significant decimal digits.
    !
    !  Machine dependent constants:
    !
    !    beta   = Radix for the floating-point system
    !    maxexp = Smallest power of beta that overflows
    !    XMAX =   Largest argument acceptable to BESI0;  Solution to
    !             equation:
    !               W(X) * (1+1/(8*X)+9/(128*X^2) = beta^maxexp
    !             where  W(X) = EXP(X)/sqrt(2*PI*X)
    !
    !    Approximate values for some important machines are:
    !
    !                             beta       maxexp       XMAX
    !
    !    CRAY-1        (S.P.)       2         8191       5682.810
    !    Cyber 180/855
    !      under NOS   (S.P.)       2         1070        745.893
    !    IEEE (IBM/XT,
    !      SUN, etc.)  (S.P.)       2          128         91.900
    !    IEEE (IBM/XT,
    !      SUN, etc.)  (D.P.)       2         1024        713.986
    !    IBM 3033      (D.P.)      16           63        178.182
    !    VAX           (S.P.)       2          127         91.203
    !    VAX D-Format  (D.P.)       2          127         91.203
    !    VAX G-Format  (D.P.)       2         1023        713.293
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 October 2008
    !
    !  Author:
    !
    !    Original FORTRAN77 version by William Cody, Laura Stoltz.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) ARG, the argument.
    !
    !    Output, real ( kind = 8 ) BESSEL_I0, the value of the modified
    !    Bessel function of the first kind.
    !
      implicit none
      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
      real ( kind = DP ) a
      real ( kind = DP ) arg
      real ( kind = DP ) b
      !real ( kind = DP ) bessel_i0
      real ( kind = DP ), parameter :: exp40 = 2.353852668370199854D+17
      integer ( kind = 4 ) i
      real ( kind = DP ), parameter, dimension ( 15 ) :: p = (/ &
        -5.2487866627945699800D-18, &
        -1.5982226675653184646D-14, &
        -2.6843448573468483278D-11, &
        -3.0517226450451067446D-08, &
        -2.5172644670688975051D-05, &
        -1.5453977791786851041D-02, &
        -7.0935347449210549190D+00, &
        -2.4125195876041896775D+03, &
        -5.9545626019847898221D+05, &
        -1.0313066708737980747D+08, &
        -1.1912746104985237192D+10, &
        -8.4925101247114157499D+11, &
        -3.2940087627407749166D+13, &
        -5.5050369673018427753D+14, &
        -2.2335582639474375249D+15 /)
      real ( kind = DP ), parameter, dimension ( 8 ) :: pp = (/ &
        -3.9843750000000000000D-01, &
         2.9205384596336793945D+00, &
        -2.4708469169133954315D+00, &
         4.7914889422856814203D-01, &
        -3.7384991926068969150D-03, &
        -2.6801520353328635310D-03, &
         9.9168777670983678974D-05, &
        -2.1877128189032726730D-06 /)
      real ( kind = DP ), parameter, dimension ( 5 ) :: q = (/ &
        -3.7277560179962773046D+03, &
         6.5158506418655165707D+06, &
        -6.5626560740833869295D+09, &
         3.7604188704092954661D+12, &
        -9.7087946179594019126D+14 /)
      real ( kind = DP ), parameter, dimension ( 7 ) :: qq = (/ &
        -3.1446690275135491500D+01, &
         8.5539563258012929600D+01, &
        -6.0228002066743340583D+01, &
         1.3982595353892851542D+01, &
        -1.1151759188741312645D+00, &
         3.2547697594819615062D-02, &
        -5.5194330231005480228D-04 /)
      real ( kind = DP ), parameter :: rec15 = 6.6666666666666666666D-02
      real ( kind = DP ) sump
      real ( kind = DP ) sumq
      real ( kind = DP ) value
      real ( kind = DP ) x
      real ( kind = DP ), parameter :: xmax = 91.9D+00
      real ( kind = DP ) xx

      x = abs ( arg )

      if ( x < epsilon ( arg ) ) then
        value = 1.0D+00
      else if ( x < 15.0D+00 ) then
    !
    !  EPSILON ( ARG ) <= ABS(ARG) < 15.0D+00
    !
        xx = x * x
        sump = p(1)
        do i = 2, 15
          sump = sump * xx + p(i)
        end do

        xx = xx - 225.0D+00
        sumq = (((( &
               xx + q(1) ) &
             * xx + q(2) ) &
             * xx + q(3) ) &
             * xx + q(4) ) &
             * xx + q(5)

        value = sump / sumq

      else if ( 15.0D+00 <= x ) then

        if ( xmax < x ) then
          value = huge ( value )
        else
    !
    !  15.0D+00 <= ABS(ARG)
    !
          xx = 1.0D+00 / x - rec15

          sump = ((((((  &
                      pp(1) &
               * xx + pp(2) ) &
               * xx + pp(3) ) &
               * xx + pp(4) ) &
               * xx + pp(5) ) &
               * xx + pp(6) ) &
               * xx + pp(7) ) &
               * xx + pp(8)

          sumq = (((((( &
                 xx + qq(1) ) &
               * xx + qq(2) ) &
               * xx + qq(3) ) &
               * xx + qq(4) ) &
               * xx + qq(5) ) &
               * xx + qq(6) ) &
               * xx + qq(7)

          value = sump / sumq
    !
    !  Calculation reformulated to avoid premature overflow.
    !
          if ( x <= xmax - 15.0D+00 ) then
            a = exp ( x )
            b = 1.0D+00
          else
            a = exp ( x - 40.0D+00 )
            b = exp40
          end if

          value = ( ( value * a - pp(1) * a ) / sqrt ( x ) ) * b

        end if

      end if

      bessel_i0 = value

      return
    endfunction bessel_i0

 
        


end subroutine network_sim 
