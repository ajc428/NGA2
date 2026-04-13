!> Models and feature detection for VOF
module detection_class
   use precision,      only: WP
   use config_class,   only: config
   use cclabel_class,  only: cclabel
   use string,         only: str_medium
   use iterator_class, only: iterator
   use tpns_class,     only: tpns
   use lpt_class,      only: lpt
   use partmesh_class, only: partmesh
   use monitor_class,  only: monitor
   use vfs_data_class, only: vfs_base
   use irl_fortran_interface
   implicit none
   private

   public :: detection

   type :: detection
      class(config), pointer :: cfg                                !< This is the config the solver is build for
      class(vfs_base), pointer :: vf                               !< This is the VOF object
      class(tpns), pointer :: fs                                   !< This is the flow solver object
      class(lpt), pointer :: lp_spray                              !< This is the particle solver object
      integer, dimension(:,:,:), allocatable :: liquid_gas_flip
      real(WP), dimension(:,:,:), allocatable :: local_thickness_recon
      integer, dimension(:,:,:), allocatable :: local_struct_type_recon
      integer, dimension(:,:,:), allocatable :: struct_type        !< Local feature type
      real(WP), dimension(:,:,:), allocatable :: film_edge_sensor  !< Edge sensing (higher is edge)
      integer, dimension(:,:,:), allocatable :: lig_edge_sensor    !< Edge sensing (higher is edge)
      real(WP), dimension(:,:,:,:), allocatable :: edge_normal     !< Edge normal
      integer, dimension(:,:,:), allocatable :: recon_type         !< Reconstruction type
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi          !< Cell-centered velocities
      real(WP), dimension(:,:), allocatable :: lig_timers
      integer, dimension(:,:,:), allocatable :: old_id
      type(cclabel)     :: ccl_drop,ccl_film,ccl_lig,ccl_thin,ccl_recon     !< CCL structures


      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      real(WP) :: vof_removed              !< Integral of VOF removed
      integer  :: nlayer=4                 !< Size of buffer layer for VOF removal

      !> Drop transfer modeling
      logical :: use_drop_transfer !< Do we use droplet transfer
      logical :: use_film_transfer !< Do we use film transfer
      logical :: use_lig_transfer  !< Do we use ligament transfer
      logical :: use_secondary     !< Do we use secondary breakup
      type(lpt)      :: lp         !< Lagrangian particle tracking
      type(monitor)  :: pfile      !< Particle monitoring
      type(partmesh) :: pmesh      !< Particle mesh for lpt
      real(WP) :: dmax             !< Maximum diameter for transfer
      real(WP) :: dmin             !< Minimum diameter below which transfer is automatic
      real(WP) :: ddel             !< Minimum diameter below which structure is directly deleted
      real(WP) :: emax             !< Maximum eccentricity for transfer
      real(WP) :: vof_tf_drop      !< Integral of VOF transfered by conversion to droplet
      real(WP) :: vof_deleted      !< Integral of VOF deleted
      integer  :: np_drop

      real(WP) :: frp
      real(WP) :: fmin
      real(WP) :: fd0
      real(WP) :: fbvol2dvol
      real(WP) :: vof_tf_film       
      integer  :: np_film

      real(WP) :: dw
      real(WP) :: size_ratio
      real(WP) :: vof_tf_lig
      integer  :: np_lig
      integer :: num_old_id = 0

   contains
      procedure :: initialize                !< Initialize detection object

      procedure :: sense_interface           !< Calculate various surface sensors
      procedure :: get_thickness             !< Calculate multiphasic structure thickness
      procedure :: detect_thin_regions       !< Detect thin regions for reconstruction
      procedure :: select_recon_type         !< Select reconstruction
      procedure :: build_recon_ccl           !< Build CCL for reconstruction
      procedure :: detect_film_recon_regions !< Detect lig regions for reconstruction 
      procedure :: detect_ligs_recon_regions !< Detect lig regions for reconstruction        
      procedure :: detect_film_edge          !< Detect edge regions
      procedure :: detect_lig_edge           !< Detect edge regions

      procedure :: prepare_transfer          !< Prepare for conversion to drops
      procedure :: attempt_transfer          !< Attempt conversion to drops
      procedure :: transfer_drops            !< Convert spheroids to drops
      procedure :: transfer_thin_features    !< Convert ligs to drops
      procedure :: secondary_break           !< Secondary breakup of drops
   end type detection

contains

   subroutine initialize(this,cfg,vf)
      implicit none
      class(detection), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      class(vfs_base), target, intent(in) :: vf
      integer :: values(8), k, nseed
      integer, dimension(:), allocatable :: seed

      call random_seed(size = k)
      allocate(seed(1:k))
      nseed = 0
      seed(:) = nseed
      call random_seed(put = seed)
      deallocate(seed)

      ! Point to objects
      this%cfg=>cfg
      this%vf=>vf
      allocate(this%liquid_gas_flip(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%liquid_gas_flip=0.0_WP
      allocate(this%local_thickness_recon(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%local_thickness_recon=0.0_WP
      allocate(this%local_struct_type_recon(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%local_struct_type_recon=0.0_WP
      allocate(this%struct_type(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%struct_type=0.0_WP
      allocate(this%film_edge_sensor(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%film_edge_sensor=0.0_WP
      allocate(this%lig_edge_sensor(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%lig_edge_sensor=0.0_WP
      allocate(this%edge_normal(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%edge_normal=0.0_WP
      allocate(this%recon_type(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%recon_type=0.0_WP
      allocate(this%Ui(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%Vi(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%Wi(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

      call this%ccl_recon%initialize(this%cfg%pgrid,name='ccl_recon')
   end subroutine initialize

   !> Reconstruction methods

   !> Compute interface sensors
   subroutine sense_interface(this)
      implicit none
      class(detection), intent(inout) :: this
      ! Update local thickness
      call this%get_thickness()
      call this%detect_thin_regions()
      ! Identify edge regions
      !call this%detect_edge_regions()
   end subroutine sense_interface

   !> Measure local thickness of multiphasic structure
   subroutine get_thickness(this)
      use vfs_data_class, only: VFlo,VFhi
      implicit none
      class(detection), intent(inout) :: this
      real(WP), dimension(:,:,:), allocatable :: tmp
      integer :: i,j,k,ii,jj,kk
      real(WP) :: lvol,gvol,area
      ! Reset thickness
      this%vf%thickness=0.0_WP
      ! First compute thickness based on current surface and volume moments (SD and VF)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond/full cells
               if (this%vf%mask(i,j,k).ne.0) cycle
               if (this%vf%VF(i,j,k).lt.VFlo.or.this%vf%VF(i,j,k).gt.VFhi) cycle
               ! Extract thickness estimate from local phasic volumes and surface area
               lvol=0.0_WP; gvol=0.0_WP; area=0.0_WP
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  lvol=lvol+(       this%vf%VF(ii,jj,kk))*this%cfg%vol(ii,jj,kk)
                  gvol=gvol+(1.0_WP-this%vf%VF(ii,jj,kk))*this%cfg%vol(ii,jj,kk)
                  area=area+        this%vf%SD(ii,jj,kk) *this%cfg%vol(ii,jj,kk)
               end do; end do; end do
               if (area.gt.0.0_WP) this%vf%thickness(i,j,k)=2.0_WP*min(lvol,gvol)/area
            end do
         end do
      end do
      call this%cfg%sync(this%vf%thickness)
      ! Filter thickness
      allocate(tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); tmp=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond/full cells
               if (this%vf%mask(i,j,k).ne.0) cycle
               if (this%vf%VF(i,j,k).lt.VFlo.or.this%vf%VF(i,j,k).gt.VFhi) cycle
               ! Surface-average thickness
               area=0.0_WP
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  area      =area      +this%vf%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)
                  tmp(i,j,k)=tmp(i,j,k)+this%vf%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)*this%vf%thickness(ii,jj,kk)
               end do; end do; end do
               if (area.gt.0.0_WP) tmp(i,j,k)=tmp(i,j,k)/area
            end do
         end do
      end do
      call this%cfg%sync(tmp)
      this%vf%thickness=tmp
      deallocate(tmp)
   end subroutine get_thickness

   !> Detect thin regions of the interface
   subroutine detect_thin_regions(this)
      use vfs_data_class, only: VFlo,VFhi
      use mathtools, only: normalize
      implicit none
      class(detection), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk,dim,dir,ni
      integer , dimension(3) :: pos
      real(WP), dimension(3) :: n1,n2,c1,c2
      real(WP), dimension(:,:,:), allocatable :: mysensor
      real(WP) :: a1,a2
      ! Default value is 0
      this%vf%thin_sensor=0.0_WP
      ! First pass to handle 2-plane cells
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Skip wall/bcond cells
               if (this%vf%mask(i,j,k).ne.0) cycle
               ! Skip full cells
               if (this%vf%VF(i,j,k).lt.VFlo.or.this%vf%VF(i,j,k).gt.VFhi) cycle
               ! Detect thin regions from local polygon data
               n1=calculateNormal(this%vf%interface_polygon(1,i,j,k))
               if (getNumberOfVertices(this%vf%interface_polygon(2,i,j,k)).gt.0) then
                  n2=calculateNormal(this%vf%interface_polygon(2,i,j,k))
                  ! Check normal orientation to identify thin regions
                  if (dot_product(n1,n2).lt.this%vf%thin_thld_dotprod) then
                     ! Check if liquid or gas
                     c1=calculateCentroid(this%vf%interface_polygon(1,i,j,k))
                     c2=calculateCentroid(this%vf%interface_polygon(2,i,j,k))
                     if (dot_product(c2-c1,n2).gt.0.0_WP) then
                        ! Thin liquid region
                        this%vf%thin_sensor(i,j,k)=1.0_WP
                     else
                        ! Thin gas region
                        this%vf%thin_sensor(i,j,k)=2.0_WP
                     end if
                  end if
               end if
            end do
         end do
      end do
      ! Second pass to extend sensor
      allocate(mysensor(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); mysensor=this%vf%thin_sensor
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond cells
               if (this%vf%mask(i,j,k).ne.0) cycle
               ! Skip full cells
               if (this%vf%VF(i,j,k).lt.VFlo.or.this%vf%VF(i,j,k).gt.VFhi) cycle
               ! Sensor is still zero, check direct neighbors
               if (this%vf%thin_sensor(i,j,k).eq.0.0_WP) then
                  do dim=1,3
                     do dir=-1,1,2
                        pos=0; pos(dim)=dir; ii=i+pos(1); jj=j+pos(2); kk=k+pos(3)
                        if (this%vf%thin_sensor(ii,jj,kk).gt.0.0_WP) mysensor(i,j,k)=this%vf%thin_sensor(ii,jj,kk)
                     end do
                  end do
               end if
            end do
         end do
      end do
      this%vf%thin_sensor=mysensor
      call this%cfg%sync(this%vf%thin_sensor)
      deallocate(mysensor)
      ! Final pass to check single-plane cells
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond cells
               if (this%vf%mask(i,j,k).ne.0) cycle
               ! Skip full cells
               if (this%vf%VF(i,j,k).lt.VFlo.or.this%vf%VF(i,j,k).gt.VFhi) cycle
               ! Sensor is still zero and 1-plane cell, check direct neighbors
               if (this%vf%thin_sensor(i,j,k).eq.0.0_WP.and.getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)).eq.1) then
                  n1=calculateNormal(this%vf%interface_polygon(1,i,j,k))
                  do dim=1,3
                     do dir=-1,1,2
                        pos=0; pos(dim)=dir; ii=i+pos(1); jj=j+pos(2); kk=k+pos(3)
                        ! Skip full cells
                        if (this%vf%VF(ii,jj,kk).lt.VFlo.or.this%vf%VF(ii,jj,kk).gt.VFhi) cycle
                        ! Check normal orientation to identify thin regions
                        ! If neighbor has two planes, then surface-average its normals and centroids
                        if (getNumberOfVertices(this%vf%interface_polygon(2,ii,jj,kk)).gt.0) then
                           a1=calculateVolume(this%vf%interface_polygon(1,ii,jj,kk))/this%cfg%meshsize(ii,jj,kk)
                           a2=calculateVolume(this%vf%interface_polygon(2,ii,jj,kk))/this%cfg%meshsize(ii,jj,kk)
                           n2=normalize(a1*calculateNormal(this%vf%interface_polygon(1,ii,jj,kk))&
                           &           +a2*calculateNormal(this%vf%interface_polygon(2,ii,jj,kk)))
                           if (dot_product(n1,n2).lt.this%vf%thin_thld_dotprod) then
                              c1=calculateCentroid(this%vf%interface_polygon(1,i ,j ,k ))
                              c2=(a1*calculateCentroid(this%vf%interface_polygon(1,ii,jj,kk))&
                              &  +a2*calculateCentroid(this%vf%interface_polygon(2,ii,jj,kk)))/(a1+a2)
                              ! Check if liquid or gas
                              if (dot_product(c2-c1,n2).gt.0.0_WP.and.dot_product(c1-c2,n1).gt.0.0_WP) then
                                 this%vf%thin_sensor(i,j,k)=1.0_WP
                              else if (dot_product(c2-c1,n2).lt.0.0_WP.and.dot_product(c1-c2,n1).lt.0.0_WP) then
                                 this%vf%thin_sensor(i,j,k)=2.0_WP
                              else
                                 this%vf%thin_sensor(i,j,k)=this%vf%thin_sensor(ii,jj,kk) ! what if it hasn't been assigned a thin sensor value yet?                          
                              end if
                           end if
                        else
                           n2=calculateNormal(this%vf%interface_polygon(1,ii,jj,kk))
                           if (dot_product(n1,n2).lt.this%vf%thin_thld_dotprod) then
                              ! Check if liquid or gas
                              c1=calculateCentroid(this%vf%interface_polygon(1,i ,j ,k ))
                              c2=calculateCentroid(this%vf%interface_polygon(1,ii,jj,kk))
                              ! Check if liquid or gas
                              if (dot_product(c2-c1,n2).gt.0.0_WP.and.dot_product(c1-c2,n1).gt.0.0_WP) then
                                 this%vf%thin_sensor(i,j,k)=1.0_WP
                              else if (dot_product(c2-c1,n2).lt.0.0_WP.and.dot_product(c1-c2,n1).lt.0.0_WP) then
                                 this%vf%thin_sensor(i,j,k)=2.0_WP
                              else
                                 this%vf%thin_sensor(i,j,k)=this%vf%thin_sensor(ii,jj,kk) ! what if it hasn't been assigned a thin sensor value yet?                          
                              end if
                           end if
                        end if
                     end do
                  end do
               end if
            end do
         end do
      end do
      call this%cfg%sync(this%vf%thin_sensor)
      ! Finally, ensure thin regions have small enough thickness
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%vf%thickness(i,j,k).gt.this%vf%thin_thld_max*this%cfg%meshsize(i,j,k)) this%vf%thin_sensor(i,j,k)=0.0_WP
            end do
         end do
      end do
   end subroutine detect_thin_regions

   !> Select reconstruction
   subroutine select_recon_type(this)
      implicit none
      class(detection), intent(inout) :: this

      !this%recon_type = 0.0_WP
      this%struct_type = 0.0_WP
      this%lig_edge_sensor=0.0_WP
      call this%build_recon_ccl()
      call this%detect_film_recon_regions()
      call this%detect_ligs_recon_regions()
      call this%detect_lig_edge()
   end subroutine

   !> Build reconstruction CCL
   subroutine build_recon_ccl(this)
      use vfs_data_class, only: VFlo,VFhi
      use mathtools, only: pi
      use mpi_f08
      use parallel,  only: MPI_REAL_WP
      implicit none
      class(detection), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk

      this%local_thickness_recon=0.0_WP
      this%local_struct_type_recon=0
   
      call get_liginfo()
      call this%ccl_recon%build(make_label,same_label)

      contains 
      ! Calculate thickness and struct_type based on moment of inertia
      subroutine get_liginfo()
         implicit none
         real(WP) :: tmpvol,tmparea,tmpvol1,tmpvol2,fluid_vol
         real(WP), dimension(1:3) :: tmpxvol, tmpL
         integer :: nneigh_moi, nneigh_thickness
         real(WP) :: x1,x2,phi
         real(WP), dimension(3)   :: d
         real(WP), dimension(3,3) :: A
         nneigh_moi=2; nneigh_thickness=3

         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (this%vf%VF(i,j,k).le.VFlo.or.this%vf%VF(i,j,k).ge.VFhi) cycle
                  ! calculate thickness
                  tmpvol=0.0_WP; tmparea=0.0_WP
                  tmpvol1=0.0_WP
                  tmpvol2=0.0_WP

                  do kk = k-nneigh_thickness,k+nneigh_thickness
                     do jj = j-nneigh_thickness,j+nneigh_thickness
                        do ii = i-nneigh_thickness,i+nneigh_thickness
                           tmpvol1 = tmpvol1 + this%vf%VF(ii,jj,kk)
                           tmpvol2 = tmpvol2 + (1.0_WP-this%vf%VF(ii,jj,kk))
                           if (this%recon_type(ii,jj,kk).eq.2.or.this%recon_type(ii,jj,kk).eq.3.or.this%recon_type(ii,jj,kk).eq.0) then
                              tmparea = tmparea + this%vf%SD(ii,jj,kk)*this%cfg%vol(i,j,k)
                           else
                              tmparea = tmparea + this%vf%SD(ii,jj,kk)*this%cfg%vol(i,j,k)*(2.0/sqrt(pi))
                           end if
                        end do
                     end do
                  end do

                  if (tmpvol1.le.tmpvol2) then
                     tmpvol = tmpvol1*this%cfg%vol(i,j,k)
                     this%liquid_gas_flip(i,j,k)=1
                  else
                    tmpvol = tmpvol2*this%cfg%vol(i,j,k)
                    this%liquid_gas_flip(i,j,k)=-1
                  end if

                  ! Calculate thickness
                  if (this%vf%VF(i,j,k).le.VFlo.and.this%liquid_gas_flip(i,j,k).eq.1) then
                     this%local_thickness_recon(i,j,k) = 0.0_WP
                  else if (this%vf%VF(i,j,k).ge.VFhi.and.this%liquid_gas_flip(i,j,k).eq.-1) then
                     this%local_thickness_recon(i,j,k) = 0.0_WP
                  else if (tmparea .gt. 0.0_WP) then    
                     this%local_thickness_recon(i,j,k) = 2.0_WP*tmpvol/(tmparea+tiny(1.0_WP))
                  else
                     this%local_thickness_recon(i,j,k) = 3.5_WP*this%cfg%min_meshsize
                  end if

                  if (this%local_thickness_recon(i,j,k).gt.VFlo) then
                     ! Calculate moi
                     tmpvol=0.0_WP; tmpxvol=0.0_WP; A=0.0_WP
                     ! First pass to accumulate volume, surface area, and position
                     do kk = k-nneigh_moi,k+nneigh_moi
                        do jj = j-nneigh_moi,j+nneigh_moi
                           do ii = i-nneigh_moi,i+nneigh_moi
                              if (this%liquid_gas_flip(i,j,k).eq.1) then
                                 tmpvol = tmpvol + this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
                                 tmpxvol = tmpxvol + this%vf%Lbary(:,ii,jj,kk)*this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
                              else
                                 tmpvol = tmpvol + (1.0_WP-this%vf%VF(ii,jj,kk))*this%cfg%vol(i,j,k)
                                 tmpxvol = tmpxvol + this%vf%Gbary(:,ii,jj,kk)*(1.0_WP-this%vf%VF(ii,jj,kk))*this%cfg%vol(i,j,k)
                              end if
                           end do
                        end do
                     end do

                     ! Second pass to accumulate moment of inertia
                     tmpxvol = tmpxvol/tmpvol
                     do kk = k-nneigh_moi,k+nneigh_moi
                        do jj = j-nneigh_moi,j+nneigh_moi
                           do ii = i-nneigh_moi,i+nneigh_moi
                              ! Location of film node
                              if (this%liquid_gas_flip(i,j,k).eq.1) then
                                 tmpL = this%vf%Lbary(:,ii,jj,kk) - tmpxvol
                                 fluid_vol = this%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
                                 A(1,1)=A(1,1)+fluid_vol*(tmpL(2)**2+tmpL(3)**2)
                                 A(2,2)=A(2,2)+fluid_vol*(tmpL(1)**2+tmpL(3)**2)
                                 A(3,3)=A(3,3)+fluid_vol*(tmpL(1)**2+tmpL(2)**2)
                                 A(1,2)=A(1,2)-fluid_vol*tmpL(1)*tmpL(2)
                                 A(1,3)=A(1,3)-fluid_vol*tmpL(1)*tmpL(3)
                                 A(2,3)=A(2,3)-fluid_vol*tmpL(2)*tmpL(3)   
                              else
                                 tmpL = this%vf%Gbary(:,ii,jj,kk) - tmpxvol
                                 fluid_vol = this%cfg%vol(ii,jj,kk)*(1.0_WP-this%vf%VF(ii,jj,kk))
                                 A(1,1)=A(1,1)+fluid_vol*(tmpL(2)**2+tmpL(3)**2)
                                 A(2,2)=A(2,2)+fluid_vol*(tmpL(1)**2+tmpL(3)**2)
                                 A(3,3)=A(3,3)+fluid_vol*(tmpL(1)**2+tmpL(2)**2)
                                 A(1,2)=A(1,2)-fluid_vol*tmpL(1)*tmpL(2)
                                 A(1,3)=A(1,3)-fluid_vol*tmpL(1)*tmpL(3)
                                 A(2,3)=A(2,3)-fluid_vol*tmpL(2)*tmpL(3)  
                              end if
                           end do
                        end do
                     end do
                     A(2,1) = A(1,2)
                     A(3,1) = A(1,3)
                     A(3,2) = A(2,3)

                     x1 = A(1,1)**2+A(2,2)**2+A(3,3)**2-A(1,1)*A(2,2)-A(1,1)*A(3,3)-A(2,2)*A(3,3)+3*(A(1,2)**2+A(1,3)**2+A(2,3)**2)
                     x2 = -(2*A(1,1)-A(2,2)-A(3,3))*(2*A(2,2)-A(1,1)-A(3,3))*(2*A(3,3)-A(1,1)-A(2,2))+9.0_WP*((2*A(3,3)-A(1,1)-A(2,2))*A(1,2)**2+(2*A(2,2)-A(1,1)-A(3,3))*A(1,3)**2+(2*A(1,1)-A(2,2)-A(3,3))*A(2,3)**2)-54.0_WP*A(1,2)*A(1,3)*A(2,3)

                     phi = atan2(sqrt(max(0.0_WP, 4*x1**3 - x2**2)), x2)

                     d(1) = (A(1,1)+A(2,2)+A(3,3)-2*sqrt(x1)*cos(phi/3.0_WP))/3.0_WP
                     d(2) = (A(1,1)+A(2,2)+A(3,3)+2*sqrt(x1)*cos((phi+pi)/3.0_WP))/3.0_WP
                     d(3) = (A(1,1)+A(2,2)+A(3,3)+2*sqrt(x1)*cos((phi-pi)/3.0_WP))/3.0_WP

                     ! Calculate local struct type
                     !call dsyev('V','U',3,A,3,d,work,lwork,info)
                     d=max(0.0_WP,d)
                     if ((d(3).lt.1.25_WP*d(2)).and.(d(2).gt.1.25_WP*d(1)) ) then
                        this%local_struct_type_recon(i,j,k) = 1
                     else if ( (d(3) .gt. 1.25_WP * d(2)) .and. (d(2) .lt. 1.25_WP * d(1)) ) then
                        this%local_struct_type_recon(i,j,k) = 2
                     end if
                  end if
               end do 
            end do 
         end do
         this%recon_type = 0.0_WP
         call this%cfg%sync(this%recon_type)
         call this%cfg%sync(this%local_thickness_recon)
         call this%cfg%sync(this%local_struct_type_recon)
         call this%cfg%sync(this%liquid_gas_flip)
      end subroutine get_liginfo
   
      !> Function that identifies cells that need a label
      logical function make_label(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         if ((this%local_thickness_recon(i,j,k).lt.1.0_WP*this%cfg%min_meshsize)) then
            make_label=.true.
         else
            make_label=.false.
         end if
      end function make_label
   
      !> Function that identifies if cell pairs have same label
      logical function same_label(i1,j1,k1,i2,j2,k2)
         implicit none
         integer, intent(in) :: i1,j1,k1,i2,j2,k2
         if (this%liquid_gas_flip(i1,j1,k1).eq.this%liquid_gas_flip(i2,j2,k2).and.this%local_struct_type_recon(i1,j1,k1).eq.this%local_struct_type_recon(i2,j2,k2)) then
            same_label=.true.
         else
            same_label=.false.
         end if
      end function same_label
   end subroutine build_recon_ccl

   !> Detect film structures for reconstruction
   subroutine detect_film_recon_regions(this)
      use vfs_data_class, only: VFlo,VFhi
      use mathtools, only: pi, normalize
      use mpi_f08
      use parallel,  only: MPI_REAL_WP
      implicit none
      class(detection), intent(inout) :: this

      integer(IRL_SignedIndex_t) :: i,j,k
      integer :: ind,ii,jj,kk
      type(VMAN_type) :: volume_moments_and_normal
      real(WP) :: surface_area,dot_result,surf_dot_pos_sum,surf_dot_neg_sum
      real(WP), dimension(3) :: surface_norm
      real(WP), dimension(:,:,:), allocatable :: tmp,tmp1
      real(WP), dimension(:,:), allocatable :: normals_adj
      real(WP), dimension(:), allocatable :: area_adj
      real(WP), dimension(:), allocatable :: norm_pos_loc,norm_neg_loc
      integer :: n,nn,size_adj,current_size,new_size
      real(WP), dimension(:,:,:), allocatable :: norm_pos
      real(WP), dimension(:,:,:), allocatable :: norm_neg

      call new(volume_moments_and_normal)
      current_size = 100
      allocate(normals_adj (1:current_size,1:3)); normals_adj =0.0_WP
      allocate(area_adj    (1:current_size));     area_adj    =0.0_WP
      allocate(norm_pos_loc(1:current_size));     norm_pos_loc=0.0_WP
      allocate(norm_neg_loc(1:current_size));     norm_neg_loc=0.0_WP

      ! Zonghao's colinearity metric
      allocate(norm_pos(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); norm_pos=0.0_WP
      allocate(norm_neg(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); norm_neg=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond/full cells
               if (this%vf%mask(i,j,k).ne.0) cycle
               if (this%vf%VF(i,j,k).lt.VFlo.or.this%vf%VF(i,j,k).gt.VFhi) cycle
               ! Count the number of triangles
               surface_area=0.0_WP; surface_norm=0.0_WP; new_size=0
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  new_size=new_size+getSize(this%vf%triangle_moments_storage(ii,jj,kk))
               end do; end do; end do
               ! Allocate local storage
               if (new_size.gt.current_size) then
                  deallocate(normals_adj, area_adj, norm_pos_loc, norm_neg_loc)
                  current_size = new_size * 1.5
                  allocate(normals_adj (1:current_size,1:3)); normals_adj =0.0_WP
                  allocate(area_adj    (1:current_size));     area_adj    =0.0_WP
                  allocate(norm_pos_loc(1:current_size));     norm_pos_loc=0.0_WP
                  allocate(norm_neg_loc(1:current_size));     norm_neg_loc=0.0_WP
               else
                  norm_pos_loc(1:new_size) = 0.0_WP
                  norm_neg_loc(1:new_size) = 0.0_WP
               end if
               ! Get surface area and normals of each triangle
               size_adj=0
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  do ind=0,getSize(this%vf%triangle_moments_storage(ii,jj,kk))-1
                     call getMoments(this%vf%triangle_moments_storage(ii,jj,kk),ind,volume_moments_and_normal)
                     size_adj=size_adj+1
                     area_adj(size_adj)     =getVolume(volume_moments_and_normal)
                     normals_adj(size_adj,:)=normalize(getNormal(volume_moments_and_normal))
                  end do
               end do; end do; end do
               surface_area=sum(area_adj(1:size_adj))
               if (surface_area.gt.0.0_WP) then
                  surf_dot_pos_sum=0.0_WP; surf_dot_neg_sum=0.0_WP
                  ! Get the postive and negative projected surface area
                  do n=1,size_adj
                     do nn=1,size_adj
                        if (n.eq.nn) cycle
                        dot_result=dot_product(normals_adj(n,:),normals_adj(nn,:))
                        if (dot_result.ge.0.0_WP) norm_pos_loc(n)=norm_pos_loc(n)+area_adj(nn)*dot_result
                        if (dot_result.lt.0.0_WP) norm_neg_loc(n)=norm_neg_loc(n)-area_adj(nn)*dot_result
                     end do
                     if ((surface_area - area_adj(n)) .gt. 1.0e-12_WP) then
                        norm_pos_loc(n)=norm_pos_loc(n)/(surface_area-area_adj(n))
                        norm_neg_loc(n)=norm_neg_loc(n)/(surface_area-area_adj(n))
                     end if
                  end do
                  ! Get the norms based on surface area weighting of the projected surface area
                  do n=1,size_adj
                     surf_dot_pos_sum=surf_dot_pos_sum+norm_pos_loc(n)*area_adj(n)
                     surf_dot_neg_sum=surf_dot_neg_sum+norm_neg_loc(n)*area_adj(n)
                  end do
                  norm_pos(i,j,k)=surf_dot_pos_sum/surface_area
                  norm_neg(i,j,k)=surf_dot_neg_sum/surface_area
               end if
            end do
         end do
      end do
      call this%cfg%sync(norm_pos);call this%cfg%sync(norm_neg)
      ! Filter metric
      allocate(tmp (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); tmp =0.0_WP
      allocate(tmp1(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); tmp1=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond/full cells
               if (this%vf%mask(i,j,k).ne.0) cycle
               if (this%vf%VF(i,j,k).lt.VFlo.or.this%vf%VF(i,j,k).gt.VFhi) cycle
               ! Surface-averaged normal magnitude
               surface_area=0.0_WP
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  surface_area=surface_area+this%vf%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)
                  tmp(i,j,k)  =tmp(i,j,k)  +this%vf%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)*norm_pos(ii,jj,kk)
                  tmp1(i,j,k) =tmp1(i,j,k) +this%vf%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)*norm_neg(ii,jj,kk)
               end do; end do; end do
               if (surface_area.gt.0.0_WP) then
                  tmp(i,j,k) =tmp(i,j,k) /surface_area
                  tmp1(i,j,k)=tmp1(i,j,k)/surface_area
               end if
               if (this%struct_type(i,j,k).ne.1 .and. .not.((tmp(i,j,k)-tmp1(i,j,k)).ge.0.5_WP.or.(((tmp(i,j,k)-tmp1(i,j,k)).lt.0.5_WP).and.(tmp(i,j,k)+tmp1(i,j,k).lt.0.75_WP)))) then
                  this%struct_type(i,j,k) = 3
                  this%recon_type(i,j,k) = 3
               else if (this%struct_type(i,j,k).ne.1) then
                  this%struct_type(i,j,k) = 2
                  this%recon_type(i,j,k) = 2
               end if
            end do
         end do
      end do
      deallocate(tmp,tmp1)
      deallocate(norm_pos,norm_neg)
      deallocate(normals_adj, area_adj, norm_pos_loc, norm_neg_loc)
   end subroutine detect_film_recon_regions

   !> Detect ligament structures for reconstruction
   subroutine detect_ligs_recon_regions(this)
      use vfs_data_class, only: VFlo,VFhi
      use mathtools, only: pi
      use mpi_f08
      use parallel,  only: MPI_REAL_WP
      implicit none
      class(detection), intent(inout) :: this
      real(WP), dimension(:)   , allocatable :: lthc
      real(WP), dimension(:)   , allocatable :: lnum
      real(WP), dimension(:)   , allocatable :: lper
      integer :: n,m,ierr,i,j,k,ii,jj,kk
      
      if (this%ccl_recon%nstruct.ge.1) then
   
         ! Allocate ligament stats arrays
         allocate(lthc(1:this%ccl_recon%nstruct)); lthc=HUGE(1.0_WP)
         allocate(lnum(1:this%ccl_recon%nstruct)); lnum=0.0_WP
         allocate(lper(1:this%ccl_recon%nstruct)); lper=0.0_WP

         ! First pass to accumulate number of cells, min thickness and ligament percentage
         do n=1,this%ccl_recon%nstruct
            ! Loop over cells in structure
            lnum(n)=lnum(n)+1.0_WP*this%ccl_recon%struct(n)%n_
            do m=1,this%ccl_recon%struct(n)%n_
               ! Get cell indices
               i=this%ccl_recon%struct(n)%map(1,m)
               j=this%ccl_recon%struct(n)%map(2,m)
               k=this%ccl_recon%struct(n)%map(3,m)

               lthc(n)=min(lthc(n),this%local_thickness_recon(i,j,k))
               !lthc(n)=lthc(n)+this%local_thickness_recon(i,j,k)
               if (this%local_struct_type_recon(i,j,k).eq.1) lper(n)=lper(n)+1.0_WP
            end do
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,lthc,1*this%ccl_recon%nstruct,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,lnum,1*this%ccl_recon%nstruct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,lper,1*this%ccl_recon%nstruct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      
         ! Second pass to set final struct_type for reconstruction
         do n=1,this%ccl_recon%nstruct
            ! Calculate the percentage of ligament structure type
            lper(n)=lper(n)/lnum(n)
            
            if ((lthc(n).le.0.8*this%cfg%min_meshsize).and.(lper(n).ge.0.9_WP).and.(lnum(n).ge.3)) then
               do m=1,this%ccl_recon%struct(n)%n_
                  i=this%ccl_recon%struct(n)%map(1,m)
                  j=this%ccl_recon%struct(n)%map(2,m)
                  k=this%ccl_recon%struct(n)%map(3,m)
                  if (this%struct_type(i,j,k).ne.3 .and. this%local_struct_type_recon(i,j,k).eq.1 .and. this%local_thickness_recon(i,j,k).le.0.8*this%cfg%min_meshsize .and. this%lig_edge_sensor(i,j,k).lt.1) then 
                     this%struct_type(i,j,k) = 1
                     this%recon_type(i,j,k) = 1
                  else if (this%struct_type(i,j,k).ne.3) then
                     this%struct_type(i,j,k) = 2
                     this%recon_type(i,j,k) = 2
                  end if
               end do
            else
               do m=1,this%ccl_recon%struct(n)%n_
                  i=this%ccl_recon%struct(n)%map(1,m)
                  j=this%ccl_recon%struct(n)%map(2,m)
                  k=this%ccl_recon%struct(n)%map(3,m)
                  if (this%struct_type(i,j,k).ne.3) then
                     this%struct_type(i,j,k) = 2
                     this%recon_type(i,j,k) = 2
                  end if
               end do
            end if
         end do
         deallocate(lthc,lnum,lper)
      end if
      call this%cfg%sync(this%struct_type)
   end subroutine detect_ligs_recon_regions

   !> Detect edge-like regions of the interface
   subroutine detect_film_edge(this)
      use vfs_data_class, only: VFlo,VFhi
      implicit none
      class(detection), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk,ni
      real(WP) :: fvol,myvol,volume_sensor,surface_sensor
      real(WP), dimension(3) :: fbary,mybary
      real(WP), dimension(:,:,:)  , allocatable :: s_tmp
      real(WP), dimension(:,:,:,:), allocatable :: v_tmp
      real(WP) :: surface_area
      ! Default value is 0
      this%film_edge_sensor=0.0_WP
      this%edge_normal=0.0_WP
      ! Traverse domain and compute sensors
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond cells
               if (this%vf%mask(i,j,k).ne.0) cycle
               ! Skip full cells
               if (this%vf%VF(i,j,k).lt.VFlo.or.this%vf%VF(i,j,k).gt.VFhi) cycle
               ! Compute filtered volume barycenter sensor
               fvol=0.0_WP; fbary=0.0_WP
               do kk=k-2,k+2
                  do jj=j-2,j+2
                     do ii=i-2,i+2
                        myvol=this%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
                        fvol =fvol +myvol
                        fbary=fbary+myvol*this%vf%Lbary(:,ii,jj,kk)
                     end do
                  end do
               end do
               fbary=fbary/fvol
               volume_sensor=norm2(fbary-this%vf%Lbary(:,i,j,k))/this%cfg%meshsize(i,j,k)
               ! Compute filtered surface barycenter sensor
               fvol=0.0_WP; mybary=0.0_WP
               do ni=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                  if (getNumberOfVertices(this%vf%interface_polygon(ni,i,j,k)).ne.0) then
                  myvol =abs(calculateVolume(this%vf%interface_polygon(ni,i,j,k)))
                  fvol  =fvol  +myvol
                  mybary=mybary+myvol*calculateCentroid(this%vf%interface_polygon(ni,i,j,k))
                  end if
               end do
               mybary=mybary/fvol
               fvol=0.0_WP; fbary=0.0_WP
               do kk=k-2,k+2
                  do jj=j-2,j+2
                     do ii=i-2,i+2
                        do ni=1,getNumberOfPlanes(this%vf%liquid_gas_interface(ii,jj,kk))
                           if (getNumberOfVertices(this%vf%interface_polygon(ni,ii,jj,kk)).ne.0) then
                              myvol=abs(calculateVolume(this%vf%interface_polygon(ni,ii,jj,kk)))
                              fvol =fvol +myvol
                              fbary=fbary+myvol*calculateCentroid(this%vf%interface_polygon(ni,ii,jj,kk))
                           end if
                        end do
                     end do
                  end do
               end do
               fbary=fbary/fvol
               surface_sensor=norm2(fbary-mybary)/this%cfg%meshsize(i,j,k)
               ! Aggregate into an edge sensor
               this%film_edge_sensor(i,j,k)=volume_sensor*surface_sensor
               ! Clip based on thickness
               if (this%vf%thickness(i,j,k).gt.this%vf%thin_thld_max*this%cfg%meshsize(i,j,k)) this%film_edge_sensor(i,j,k)=0.0_WP
               ! Finally, store edge orientation
               if (this%film_edge_sensor(i,j,k).gt.0.0_WP) then
                  this%edge_normal(:,i,j,k)=(fbary-mybary)/(norm2(fbary-mybary)+epsilon(1.0_WP))
               end if
            end do
         end do
      end do
      ! Communicate
      call this%cfg%sync(this%film_edge_sensor)
      call this%cfg%sync(this%edge_normal)
      ! Apply an extra step of surface smoothing to our edge info
      allocate(s_tmp(    this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); s_tmp=0.0_WP
      allocate(v_tmp(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); v_tmp=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond/full cells
               if (this%vf%mask(i,j,k).ne.0) cycle
               if (this%vf%VF(i,j,k).lt.VFlo.or.this%vf%VF(i,j,k).gt.VFhi) cycle
               ! Surface-averaged normal magnitude
               surface_area=0.0_WP
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  surface_area  =  surface_area+this%vf%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)
                  s_tmp  (i,j,k)=s_tmp  (i,j,k)+this%vf%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)*this%film_edge_sensor  (ii,jj,kk)
                  v_tmp(:,i,j,k)=v_tmp(:,i,j,k)+this%vf%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)*this%edge_normal(:,ii,jj,kk)
               end do; end do; end do
               if (surface_area.gt.0.0_WP) then
                  s_tmp  (i,j,k)=s_tmp  (i,j,k)/surface_area
                  v_tmp(:,i,j,k)=v_tmp(:,i,j,k)/surface_area
                  v_tmp(:,i,j,k)=v_tmp(:,i,j,k)/(norm2(v_tmp(:,i,j,k))+epsilon(1.0_WP))
               end if
            end do
         end do
      end do
      call this%cfg%sync(s_tmp); this%film_edge_sensor=s_tmp; deallocate(s_tmp)
      call this%cfg%sync(v_tmp); this%edge_normal=v_tmp; deallocate(v_tmp)
   end subroutine detect_film_edge

   !> Detect ligament edge-like regions of the interface
   subroutine detect_lig_edge(this)
      use vfs_data_class, only: VFlo,VFhi
      implicit none
      class(detection), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk,iii,jjj,kkk,lim,n,m,di,dj,dk
      real(WP) :: count,count1
      real(IRL_double), dimension(1:3,1:3,1:3) :: counts
      real(IRL_double), dimension(0:6,-1:1,-1:1,-1:1) :: moments
      real(IRL_double) :: m000,m100,m010,m001,m000g,m100g,m010g,m001g
      real(IRL_double), dimension(0:2) :: center,centerg,diff

      ! Default value is 0
      this%lig_edge_sensor=0.0_WP
      lim=3

      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            main_loop: do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond cells
               if (this%vf%mask(i,j,k).ne.0) cycle
               ! Skip full cells
               if (this%vf%VF(i,j,k).lt.VFlo.or.this%vf%VF(i,j,k).gt.VFhi) cycle
               count = 0
               count1 = 0
               counts = 0.0_WP

               moments = 0.0_WP
               m000=0; m100=0; m010=0; m001=0
               m000g=0; m100g=0; m010g=0; m001g=0
               center=0.0_WP
               centerg=0.0_WP
               
               do kk=k-1,k+1
                  do jj=j-1,j+1
                     do ii=i-1,i+1
                        if (this%vf%mask(ii,jj,kk).ne.0.or.this%lig_edge_sensor(i,j,k).eq.2.0_WP) then 
                           this%lig_edge_sensor(i,j,k)=2.0_WP
                           cycle main_loop
                        end if
                        di = ii - i
                        dj = jj - j
                        dk = kk - k
                        if (this%ccl_recon%id(ii,jj,kk).ne.this%ccl_recon%id(i,j,k).or.(ii.eq.i.and.jj.eq.j.and.kk.eq.k)) then
                           moments(0,di,dj,dk)=0.0_WP
                           moments(1,di,dj,dk)=0.0_WP
                           moments(2,di,dj,dk)=0.0_WP
                           moments(3,di,dj,dk)=0.0_WP
                           moments(4,di,dj,dk)=0.0_WP
                           moments(5,di,dj,dk)=0.0_WP
                           moments(6,di,dj,dk)=0.0_WP
                        else
                           moments(0,di,dj,dk)=this%vf%VF(ii,jj,kk)
                           moments(1,di,dj,dk)=(this%vf%Lbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                           moments(2,di,dj,dk)=(this%vf%Lbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                           moments(3,di,dj,dk)=(this%vf%Lbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                           moments(4,di,dj,dk)=(this%vf%Gbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                           moments(5,di,dj,dk)=(this%vf%Gbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                           moments(6,di,dj,dk)=(this%vf%Gbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)

                           ! Calculate geometric moments of neighborhood
                           m000=m000+moments(0,di,dj,dk)
                           m100=m100+(moments(1,di,dj,dk)+(ii-i))*moments(0,di,dj,dk)
                           m010=m010+(moments(2,di,dj,dk)+(jj-j))*moments(0,di,dj,dk)
                           m001=m001+(moments(3,di,dj,dk)+(kk-k))*moments(0,di,dj,dk)
                           m000g=m000g+(1.0_WP-moments(0,di,dj,dk))
                           m100g=m100g+(moments(4,di,dj,dk)+(ii-i))*(1.0_WP-moments(0,di,dj,dk))
                           m010g=m010g+(moments(5,di,dj,dk)+(jj-j))*(1.0_WP-moments(0,di,dj,dk))
                           m001g=m001g+(moments(6,di,dj,dk)+(kk-k))*(1.0_WP-moments(0,di,dj,dk))
                        end if
                     end do
                  end do
               end do
               ! Calculate geometric center of neighborhood
               if (m000.gt.VFlo) center=[m100,m010,m001]/m000
               if (m000g.gt.VFlo) centerg=[m100g,m010g,m001g]/m000g
               count = 0
               count1 = 0
               
               do kk=k-1,k+1
                  do jj=j-1,j+1
                     do ii=i-1,i+1
                        di = ii - i
                        dj = jj - j
                        dk = kk - k
                        if (this%vf%mask(ii,jj,kk).ne.0.or.this%lig_edge_sensor(i,j,k).eq.2) then 
                           this%lig_edge_sensor(i,j,k)=2.0_WP
                           cycle main_loop
                        end if
                        if (this%ccl_recon%id(ii,jj,kk).ne.this%ccl_recon%id(i,j,k).or.(ii.eq.i.and.jj.eq.j.and.kk.eq.k)) then
                           count = count + 0.0_WP
                           count1 = count1 + 0.0_WP
                        else
                           count = count + ((moments(1,di,dj,dk)+di - center(0))**2)*(moments(0,di,dj,dk)) + ((moments(2,di,dj,dk)+dj - center(1))**2)*(moments(0,di,dj,dk)) + ((moments(3,di,dj,dk)+dk - center(2))**2)*(moments(0,di,dj,dk))
                           count1 = count1 + ((moments(4,di,dj,dk)+di - centerg(0))**2)*(1.0_WP-moments(0,di,dj,dk)) + ((moments(5,di,dj,dk)+dj - centerg(1))**2)*(1.0_WP-moments(0,di,dj,dk)) + ((moments(6,di,dj,dk)+dk - centerg(2))**2)*(1.0_WP-moments(0,di,dj,dk))
                        end if
                     end do
                  end do
               end do
               if (m000.gt.VFlo) count = count / m000
               if (m000g.gt.VFlo) count1 = count1 / m000g
               if (this%liquid_gas_flip(i,j,k).eq.-1) count = count1
               if (this%lig_edge_sensor(i,j,k).eq.2) then
                  this%lig_edge_sensor(i,j,k)=0.0_WP
               else if (this%struct_type(i,j,k).ne.3.and.count.le.0.2) then!.and.this%struct_type(i,j,k).eq.1.0_WP) then
                  this%lig_edge_sensor(i,j,k)=1.0_WP
                  this%recon_type(i,j,k)=2
               end if
            end do main_loop
         end do
      end do
      ! Communicate
      call this%cfg%sync(this%recon_type)
      call this%cfg%sync(this%lig_edge_sensor)

   end subroutine detect_lig_edge














   !> Breakup modeling methods

   !> Prepare transfer methods
   subroutine prepare_transfer(this,use_drop_transfer,use_film_transfer,use_lig_transfer,use_secondary,fs,lp_spray)
      use messager,  only: die
      use filesys,  only: makedir,isdir
      class(detection), intent(inout) :: this
      class(tpns), target, intent(in) :: fs
      class(lpt), target, intent(in) :: lp_spray
      logical, intent(in) :: use_drop_transfer, use_film_transfer, use_lig_transfer, use_secondary
      integer :: ierr,iunit
      character(len=str_medium) :: filename

      ! Point to objects
      this%fs=>fs
      this%lp_spray=>lp_spray

      allocate(this%lig_timers(1,1)); this%lig_timers = 0.0_WP
      allocate(this%old_id(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%old_id=0.0_WP

      this%use_drop_transfer = use_drop_transfer
      this%use_film_transfer = use_film_transfer
      this%use_lig_transfer = use_lig_transfer
      this%use_secondary = use_secondary

      ! Only initialize transfer model if used
      if (this%use_drop_transfer) then
         ! Create CCL
         call this%ccl_drop%initialize(pg=this%cfg%pgrid,name='ccl_drop')
         ! Set parameters for transfer
         this%ddel=0.2_WP*this%cfg%min_meshsize
         this%dmin=1.5_WP*this%cfg%min_meshsize
         this%dmax=10*this%cfg%min_meshsize
         this%emax=0.75_WP
         ! Zero out monitoring variables
         this%vof_tf_drop=0.0_WP
         this%np_drop=0
      end if

      if (this%use_film_transfer) then
         ! Create CCL
         call this%ccl_thin%initialize(pg=this%cfg%pgrid,name='ccl_thin')
         ! this%fmin=2.2e-6 ! Emperical minimum bag thickness from Jackiw and Ashgriz 2022
         !this%fmin=1.0e-6 ! Emperical minimum bag thickness from Jackiw and Ashgriz 2022
         this%fmin=1.0e-3
         this%fbvol2dvol=0.25_WP
         ! Zero out monitoring variables
         this%vof_tf_film=0.0_WP
         this%np_film=0
      end if

      if (this%use_lig_transfer) then
         ! Create CCL LIG
         if (.not.use_film_transfer) call this%ccl_thin%initialize(pg=this%cfg%pgrid,name='ccl_thin')
         this%dw =0.697_WP
         this%size_ratio=0.5_WP !0.015_WP!0.707_WP 
         ! Zero out monitoring variables
         this%vof_tf_lig=0.0_WP
         this%np_lig=0
         this%fmin=1.0e-3
         this%fbvol2dvol=0.25_WP
      end if

      if (this%cfg%amroot) then
         if (.not.isdir('spray-all')) call makedir('spray-all')
         filename='spray-all/droplets'
         open(newunit=iunit,file=trim(filename),form='formatted',status='unknown',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
         close(iunit)         
      end if
   end subroutine prepare_transfer

   subroutine attempt_transfer(this,dt)
      class(detection), intent(inout) :: this
      real(WP), intent(in) :: dt
      ! Zero out monitoring variables
      this%lp_spray%np_new=0
      this%lp_spray%vp_new=0.0_WP

      this%lp_spray%np_new=0
      this%lp_spray%vp_new=0.0_WP
      if (this%use_drop_transfer) call this%transfer_drops()
      if (this%use_film_transfer.or.this%use_lig_transfer) call this%transfer_thin_features(dt)
      if (this%use_secondary) call this%secondary_break()
   end subroutine attempt_transfer

   !> Transfer droplet to Lagrangian representation
   subroutine transfer_drops(this)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: pi
      use messager,  only: die
      use vfs_data_class, only: VFlo,VFhi
      class(detection), intent(inout) :: this
      real(WP), dimension(:)    , allocatable :: dvol
      real(WP), dimension(:,:)  , allocatable :: dpos
      real(WP), dimension(:,:)  , allocatable :: dvel
      real(WP), dimension(:,:,:), allocatable :: dmoi
      integer :: n,m,ierr,i,j,k,iunit,np_start
      real(WP) :: x,y,z,x0,y0,z0,diam,ecc,lmax,lmid,lmin
      character(len=str_medium) :: filename
      logical :: transfer
      ! Moment of inertia calculation using lapack
      real(WP), dimension(:), allocatable, save :: work !< Saved!
      integer, save :: lwork                            !< Saved!
      real(WP), dimension(1) :: lwork_query
      real(WP), dimension(3) :: d
      real(WP), dimension(3,3) :: A
      integer :: info
      
      ! Query optimal work array size
      if (.not.allocated(work)) then
         call dsyev('V','U',3,A,3,d,lwork_query,-1,info)
         lwork=int(lwork_query(1)); allocate(work(lwork))
      end if
      
      ! Start by performing a CCL
      call this%ccl_drop%build(make_label,same_label)
      
      ! Allocate droplet stats arrays
      allocate(dvol(1:this%ccl_drop%nstruct        )); dvol=0.0_WP
      allocate(dpos(1:this%ccl_drop%nstruct,1:3    )); dpos=0.0_WP
      allocate(dvel(1:this%ccl_drop%nstruct,1:3    )); dvel=0.0_WP
      allocate(dmoi(1:this%ccl_drop%nstruct,1:3,1:3)); dmoi=0.0_WP

      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      
      ! First pass to accumulate volume, position, and velocity
      do n=1,this%ccl_drop%nstruct
         ! Loop over cells in structure
         do m=1,this%ccl_drop%struct(n)%n_
            ! Get cell indices
            i=this%ccl_drop%struct(n)%map(1,m)
            j=this%ccl_drop%struct(n)%map(2,m)
            k=this%ccl_drop%struct(n)%map(3,m)
            ! Get cell position, accounting for periodicity
            x=this%cfg%xm(i)-this%ccl_drop%struct(n)%per(1)*this%cfg%xL
            y=this%cfg%ym(j)-this%ccl_drop%struct(n)%per(2)*this%cfg%yL
            z=this%cfg%zm(k)-this%ccl_drop%struct(n)%per(3)*this%cfg%zL
            ! Accumulate volume, position, and velocity
            dvol(n  )=dvol(n  )+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            dpos(n,:)=dpos(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[x,y,z]
            dvel(n,:)=dvel(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[this%Ui(i,j,k),this%Vi(i,j,k),this%Wi(i,j,k)]
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,dvol,1*this%ccl_drop%nstruct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,dpos,3*this%ccl_drop%nstruct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,dvel,3*this%ccl_drop%nstruct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      
      ! Second pass to accumulate moment of inertia
      do n=1,this%ccl_drop%nstruct
         ! Get drop barycenter
         x0=dpos(n,1)/dvol(n)
         y0=dpos(n,2)/dvol(n)
         z0=dpos(n,3)/dvol(n)
         ! Loop over cells in structure
         do m=1,this%ccl_drop%struct(n)%n_
            ! Get cell indices
            i=this%ccl_drop%struct(n)%map(1,m)
            j=this%ccl_drop%struct(n)%map(2,m)
            k=this%ccl_drop%struct(n)%map(3,m)
            ! Get cell position relative to drop barycenter, accounting for periodicity
            x=this%cfg%xm(i)-this%ccl_drop%struct(n)%per(1)*this%cfg%xL-x0
            y=this%cfg%ym(j)-this%ccl_drop%struct(n)%per(2)*this%cfg%yL-y0
            z=this%cfg%zm(k)-this%ccl_drop%struct(n)%per(3)*this%cfg%zL-z0
            ! Accumulate moment of inertia
            dmoi(n,1,1)=dmoi(n,1,1)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y**2+z**2)
            dmoi(n,2,2)=dmoi(n,2,2)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(z**2+x**2)
            dmoi(n,3,3)=dmoi(n,3,3)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x**2+y**2)
            dmoi(n,1,2)=dmoi(n,1,2)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*y)
            dmoi(n,1,3)=dmoi(n,1,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*z)
            dmoi(n,2,3)=dmoi(n,2,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y*z)
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,dmoi,9*this%ccl_drop%nstruct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      
      ! Third pass to generate normalized drop stats
      do n=1,this%ccl_drop%nstruct
         ! Get drop barycenter, accounting for periodicity
         dpos(n,:)=dpos(n,:)/dvol(n)
         if (this%cfg%xper.and.dpos(n,1).lt.this%cfg%x(this%cfg%imin)) dpos(n,1)=dpos(n,1)+this%cfg%xL
         if (this%cfg%yper.and.dpos(n,2).lt.this%cfg%y(this%cfg%jmin)) dpos(n,2)=dpos(n,2)+this%cfg%yL
         if (this%cfg%zper.and.dpos(n,3).lt.this%cfg%z(this%cfg%kmin)) dpos(n,3)=dpos(n,3)+this%cfg%zL
         ! Get drop velocity
         dvel(n,:)=dvel(n,:)/dvol(n)
      end do
      
      ! Zero out monitoring variables
      this%vof_tf_drop=0.0_WP
      this%vof_deleted=0.0_WP
      this%np_drop=0

      if (this%cfg%amRoot) then
         filename='spray-all/droplets'
         open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
         if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
      end if
      ! Transfer drops based on our criteria
      do n=1,this%ccl_drop%nstruct
         
         ! Compute diameter
         diam=(6.0_WP*dvol(n)/pi)**(1.0_WP/3.0_WP)
         
         ! Decide whether to transfer based on diameter
         if (diam.gt.this%dmax) then
            transfer=.false.
         else if (diam.gt.this%ddel.and.diam.le.this%dmin) then
            ! Small enough to transfer automatically
            transfer=.true.
         else
            ! In between, check eccentricity from moment of inertia tensor
            A=dmoi(n,:,:)
            call dsyev('V','U',3,A,3,d,work,lwork,info) !< On exit, A contains eigenvectors and d contains eigenvalues in ascending order
            d=max(0.0_WP,d)                             !< Get rid of very small negative values (due to machine accuracy)
            ! Get characteristic lengths of drop
            lmax=sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/dvol(n))
            lmid=sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/dvol(n))
            lmin=sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/dvol(n))
            if (lmin.eq.0.0_WP) lmin=lmid ! Handle 2D case
            ecc=sqrt(1.0_WP-lmin**2/(lmax**2+epsilon(1.0_WP)))
            if (d(3).lt.2.0*d(1)) then
               transfer=.true.
            else
               transfer=.false.
            end if
         end if
         
         ! Perform transfer
         if (transfer) then
            ! Root creates a new Lagrangian drop
            if (this%cfg%amRoot) then
               print *, "This is a drop with diam of ", diam
               np_start=this%lp_spray%np_
               ! Increment particle counter
               this%lp_spray%np_=this%lp_spray%np_+1
               ! Make room for new drop
               call this%lp_spray%resize(this%lp_spray%np_)
               ! Add the drop

               this%lp_spray%p(this%lp_spray%np_)%id  =int(1,8)
               this%lp_spray%p(this%lp_spray%np_)%d   =diam
               this%lp_spray%p(this%lp_spray%np_)%pos =dpos(n,:)
               this%lp_spray%p(this%lp_spray%np_)%vel =dvel(n,:)
               this%lp_spray%p(this%lp_spray%np_)%ind =this%cfg%get_ijk_global(dpos(n,:),[this%lp_spray%cfg%imin,this%lp_spray%cfg%jmin,this%lp_spray%cfg%kmin])
               this%lp_spray%p(this%lp_spray%np_)%flag=0
               this%lp_spray%p(this%lp_spray%np_)%dt  =0.0_WP
               this%lp_spray%p(this%lp_spray%np_)%Acol=0.0_WP
               this%lp_spray%p(this%lp_spray%np_)%Tcol=0.0_WP
               this%lp_spray%p(this%lp_spray%np_)%t   =0.0_WP

               !!! Write to droplet list !!!
               ! Output diameter, velocity, and position
               write(iunit,*) this%lp_spray%p(this%lp_spray%np_)%d,this%lp_spray%p(this%lp_spray%np_)%vel(1),this%lp_spray%p(this%lp_spray%np_)%vel(2),this%lp_spray%p(this%lp_spray%np_)%vel(3),&
               &norm2([this%lp_spray%p(this%lp_spray%np_)%vel(1),this%lp_spray%p(this%lp_spray%np_)%vel(2),this%lp_spray%p(this%lp_spray%np_)%vel(3)]),this%lp_spray%p(this%lp_spray%np_)%pos(1),&
               &this%lp_spray%p(this%lp_spray%np_)%pos(2),this%lp_spray%p(this%lp_spray%np_)%pos(3),this%lp_spray%p(this%lp_spray%np_)%id, 0  
            end if
            
            ! Zero out VF in the structure
            do m=1,this%ccl_drop%struct(n)%n_
               this%vf%VF(this%ccl_drop%struct(n)%map(1,m),this%ccl_drop%struct(n)%map(2,m),this%ccl_drop%struct(n)%map(3,m))=0.0_WP
            end do
            
            ! Increment monitoring variables
            this%vof_tf_drop=this%vof_tf_drop+dvol(n)
            this%np_drop=this%np_drop+1
            this%lp_spray%np_new=this%lp_spray%np_new+1
            this%lp_spray%vp_new=this%lp_spray%vp_new+dvol(n)
         end if
      end do

      if (this%cfg%amRoot) close(iunit)
      
      ! Synchronize VF fields
      call this%vf%sync_interface()
      call this%vf%clean_irl_and_band()
      
      ! Synchronize particles
      call this%lp_spray%sync()
      
      ! Deallocate all but work array
      deallocate(dvol,dpos,dvel,dmoi)
      
   contains
      !> Function that identifies cells that need a label
      logical function make_label(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         if (this%vf%VF(i,j,k).gt.VFlo) then
         make_label=.true.
         else
         make_label=.false.
         end if
      end function make_label
      
      !> Function that identifies if cell pairs have same label
      logical function same_label(i1,j1,k1,i2,j2,k2)
         implicit none
         integer, intent(in) :: i1,j1,k1,i2,j2,k2
         same_label=.true.
      end function same_label
      
   end subroutine transfer_drops

   subroutine transfer_thin_features(this,dt)
      use vfs_data_class, only: VFlo,VFhi
      use mathtools, only: pi,twoPi,normalize,cross_product
      use parallel,  only: MPI_REAL_WP
      use messager,  only: die
      use random,    only: random_uniform,random_gamma
      use mpi_f08
      use irl_fortran_interface
      implicit none
      class(detection), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(:)   , allocatable :: lvol
      real(WP), dimension(:)   , allocatable :: lthc
      real(WP), dimension(:)   , allocatable :: llen
      real(WP), dimension(:)   , allocatable :: lnum
      real(WP), dimension(:)   , allocatable :: lper
      real(WP), dimension(:,:) , allocatable :: lpos
      real(WP), dimension(:,:) , allocatable :: lvel
      real(WP), dimension(:,:,:), allocatable :: lmoi
      real(WP), dimension(:)   , allocatable :: lSR
      real(WP), dimension(:)   , allocatable :: xmin,xmax,ymin,ymax,zmin,zmax
      integer :: n,m,ierr,i,j,k,l,ii,jj,kk,iunit,ncell_,tmp_id,totalnewp,np_start,np_old,ip,rank
      real(WP) :: x,y,z,x0,y0,z0,lmax,lmid,lmin,d_0,d_avg,nu,d_break
      real(WP)  :: tmp_ke,curv_sum,ncurv,Vt,Vl,Vd,alpha,beta,fd0,frp
      real(WP), dimension(3) :: nref,tref,sref
      logical :: sampled, frem_active, at_edge, at_edge_total
   
      character(len=str_medium) :: filename
      real(WP), dimension(:), allocatable :: sort_ke
      integer, dimension(:), allocatable ::  sort_id,plist,dispels,lig_film
      real(WP), dimension(:,:), allocatable :: pinfo,pinfo_
      real(WP) :: minor_radius,diam,Vrim,Lrim
      real(WP) :: Trp,Lrp,Tsr,SR_tmp,T_count
      real(WP), dimension(1:3) :: tangent
      real(WP), dimension(:,:,:,:), allocatable :: SR
      integer  :: nmain,nsat
      real(WP), dimension(:,:,:), allocatable :: local_thickness
      integer,  dimension(:,:,:), allocatable :: local_struct_type
      ! Moment of inertia calculation using lapack
      real(WP), dimension(:), allocatable, save :: work !< Saved!
      integer, save :: lwork                      !< Saved!
      real(WP), dimension(1) :: lwork_query
      real(WP), dimension(3) :: d
      real(WP), dimension(3,3) :: A
      integer :: info
      real(IRL_double), dimension(1:3) :: a_aligned_Cylinder
   
      integer :: count
      real(WP), dimension(:,:), allocatable :: points
      real(WP) :: prob
      real(WP), dimension(:), allocatable :: lig_sizes
   
      real(WP), dimension(:,:), allocatable :: lig_timers_new
      integer, dimension(:,:), allocatable :: overlap
      integer, dimension(:), allocatable :: survivor_map
      integer, dimension(:,:), allocatable :: packed_lig_timers
      integer :: survivor_count, current_id
   
      type :: spline_info
         integer :: n_knots_x, n_knots_y, n_knots_z, flag
         real(WP) :: length
         real(WP), dimension(:), allocatable :: t_knots_x, c_coeffs_x
         real(WP), dimension(:), allocatable :: t_knots_y, c_coeffs_y
         real(WP), dimension(:), allocatable :: t_knots_z, c_coeffs_z
      end type spline_info
   
      type(spline_info) :: s_info
   
      ! Start by performing a CCL based on ligament criteria
      ! Get thickness and local struct_type for global information calculation
      allocate(local_thickness(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));local_thickness=0.0_WP
      allocate(local_struct_type(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));local_struct_type=0
      call get_structinfo()
      call this%ccl_thin%build(make_label,same_label)
   
      if (this%ccl_thin%nstruct.ge.1) then

         ! Record initial droplets in each processor for future outputing purpose
         np_start=this%lp_spray%np_; sampled=.false.
   
         ! Query optimal work array size
         if (.not.allocated(work)) then
            call dsyev('V','U',3,A,3,d,lwork_query,-1,info)
            lwork=int(lwork_query(1)); allocate(work(lwork))
         end if
   
         ! Allocate ligament stats arrays
         allocate(lig_film(1:this%ccl_thin%nstruct)); lig_film=0
         allocate(lvol(1:this%ccl_thin%nstruct)); lvol=0.0_WP
         allocate(lthc(1:this%ccl_thin%nstruct)); lthc=HUGE(1.0_WP)
         allocate(llen(1:this%ccl_thin%nstruct)); llen=0.0_WP
         allocate(lnum(1:this%ccl_thin%nstruct)); lnum=0.0_WP
         allocate(lper(1:this%ccl_thin%nstruct)); lper=0.0_WP
         allocate(lpos(1:this%ccl_thin%nstruct,1:3)); lpos=0.0_WP
         allocate(lvel(1:this%ccl_thin%nstruct,1:3)); lvel=0.0_WP
         allocate(lmoi(1:this%ccl_thin%nstruct,1:3,1:3)); lmoi=0.0_WP
         allocate(lSR(1:this%ccl_thin%nstruct)); lSR=-HUGE(1.0_WP)
         allocate(xmin(1:this%ccl_thin%nstruct),xmax(1:this%ccl_thin%nstruct)); xmin=HUGE(1.0_WP);xmax=-HUGE(1.0_WP)
         allocate(ymin(1:this%ccl_thin%nstruct),ymax(1:this%ccl_thin%nstruct)); ymin=HUGE(1.0_WP);ymax=-HUGE(1.0_WP)
         allocate(zmin(1:this%ccl_thin%nstruct),zmax(1:this%ccl_thin%nstruct)); zmin=HUGE(1.0_WP);zmax=-HUGE(1.0_WP)
         allocate(SR(1:6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));SR=0.0_WP
         allocate(lig_timers_new(this%ccl_thin%nstruct,2));lig_timers_new=0.0_WP
         call this%fs%get_strainrate(SR)
   
         call this%cfg%sync(this%vf%SD)
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
   
         ! First pass to accumulate volume, position, min thickness and ligament percentage
         do n=1,this%ccl_thin%nstruct
            ! Loop over cells in structure
            lnum(n)=lnum(n)+1.0_WP*this%ccl_thin%struct(n)%n_
            do m=1,this%ccl_thin%struct(n)%n_
               ! Get cell indices
               i=this%ccl_thin%struct(n)%map(1,m)
               j=this%ccl_thin%struct(n)%map(2,m)
               k=this%ccl_thin%struct(n)%map(3,m)
               ! Get cell position, accounting for periodicity
               x=this%cfg%xm(i)-this%ccl_thin%struct(n)%per(1)*this%cfg%xL
               y=this%cfg%ym(j)-this%ccl_thin%struct(n)%per(2)*this%cfg%yL
               z=this%cfg%zm(k)-this%ccl_thin%struct(n)%per(3)*this%cfg%zL
   
               ! Accumulate volume, position, velocity. Get min thickness and ligament percentage
               lvol(n)=lvol(n)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
               lpos(n,:)=lpos(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[x,y,z]
               lvel(n,:)=lvel(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[this%Ui(i,j,k),this%Vi(i,j,k),this%Wi(i,j,k)]
               lthc(n)=min(lthc(n),local_thickness(i,j,k))
               if (local_struct_type(i,j,k).eq.1) lper(n)=lper(n)+1.0_WP
   
               ! Get the structure's bounding box
               if (this%recon_type(i,j,k).eq.2) then
                  do l=1,2
                     if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
                        d = calculateCentroid(this%vf%interface_polygon(l,i,j,k))
                        xmin(n)=min(xmin(n),d(1)); xmax(n)=max(xmax(n),d(1))
                        ymin(n)=min(ymin(n),d(2)); ymax(n)=max(ymax(n),d(2))
                        zmin(n)=min(zmin(n),d(3)); zmax(n)=max(zmax(n),d(3))
                     end if
                  end do
               else if (this%recon_type(i,j,k).eq.1) then
                  if (this%vf%VF(i,j,k).ge.VFlo) then
                     d = this%vf%LBARY(:,i,j,k)
                     xmin(n)=min(xmin(n),d(1)); xmax(n)=max(xmax(n),d(1))
                     ymin(n)=min(ymin(n),d(2)); ymax(n)=max(ymax(n),d(2))
                     zmin(n)=min(zmin(n),d(3)); zmax(n)=max(zmax(n),d(3))
                  end if
               end if
            end do
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,lvol,1*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,lpos,3*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,lvel,3*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,lthc,1*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,lnum,1*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,lper,1*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,xmin,1*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,ymin,1*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,zmin,1*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,xmax,1*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,ymax,1*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,zmax,1*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
         
         ! Zero out monitoring variables
         this%vof_tf_film=0.0_WP
         this%np_film=0

         ! Second pass to accumulate moment of inertia
         do n=1,this%ccl_thin%nstruct
            ! Get ligament barycenter
            if (lvol(n) .le. VFlo) cycle
            x0=lpos(n,1)/lvol(n)
            y0=lpos(n,2)/lvol(n)
            z0=lpos(n,3)/lvol(n)
            ! Loop over cells in structure
            do m=1,this%ccl_thin%struct(n)%n_
               ! Get cell indices
               i=this%ccl_thin%struct(n)%map(1,m)
               j=this%ccl_thin%struct(n)%map(2,m)
               k=this%ccl_thin%struct(n)%map(3,m)
               ! Get cell position relative to drop barycenter, accounting for periodicity
               x=this%cfg%xm(i)-x0
               y=this%cfg%ym(j)-y0
               z=this%cfg%zm(k)-z0
               x=x-this%ccl_thin%struct(n)%per(1)*this%cfg%xL
               y=y-this%ccl_thin%struct(n)%per(2)*this%cfg%yL
               z=z-this%ccl_thin%struct(n)%per(3)*this%cfg%zL
   
               ! Accumulate moment of inertia
               lmoi(n,1,1)=lmoi(n,1,1)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y**2+z**2)
               lmoi(n,2,2)=lmoi(n,2,2)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(z**2+x**2)
               lmoi(n,3,3)=lmoi(n,3,3)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x**2+y**2)
               lmoi(n,1,2)=lmoi(n,1,2)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*y)
               lmoi(n,1,3)=lmoi(n,1,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*z)
               lmoi(n,2,3)=lmoi(n,2,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y*z)
            end do
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,lmoi,9*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      
         ! Third pass to generalize ligament stats
         do n=1,this%ccl_thin%nstruct
            if (lvol(n) .le. VFlo) cycle
            ! Get ligament barycenter, accounting for periodicity
            lpos(n,:)=lpos(n,:)/lvol(n)
            if (this%cfg%xper.and.lpos(n,1).lt.this%cfg%x(this%cfg%imin)) lpos(n,1)=lpos(n,1)+this%cfg%xL
            if (this%cfg%yper.and.lpos(n,2).lt.this%cfg%y(this%cfg%jmin)) lpos(n,2)=lpos(n,2)+this%cfg%yL
            if (this%cfg%zper.and.lpos(n,3).lt.this%cfg%z(this%cfg%kmin)) lpos(n,3)=lpos(n,3)+this%cfg%zL
            ! Get drop velocity
            lvel(n,:)=lvel(n,:)/lvol(n)
            ! Calculate the percentage of ligament structure type
            lper(n)=lper(n)/lnum(n)
            ! Calculate maximum length of the structure
            A=lmoi(n,:,:)
            call dsyev('V','U',3,A,3,d,work,lwork,info) !< On exit, A contains eigenvectors and d contains eigenvalues in ascending order
            d=max(0.0_WP,d)
            ! Replace with corrected eigenvectors for future ligament droplet placement
            lmoi(n,:,:)=A
            ! Get characteristic lengths of drop
            lmax=sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/lvol(n))
            lmid=sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/lvol(n))
            lmin=sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/lvol(n))
            if (lmin.eq.0.0_WP) lmin=lmid ! Handle 2D case
            ! Use max of bounding box and MoI-derived lengths as length
            llen(n) = max(sqrt((xmax(n)-xmin(n))**2+(ymax(n)-ymin(n))**2+(zmax(n)-zmin(n))**2),lmax)
      
            ! With the tangent direction of the ligament, we can evaluate the strain rate of each cell of the ligament
            tangent = lmoi(n,:,1)
            do m=1,this%ccl_thin%struct(n)%n_
               ! Get cell indices
               i=this%ccl_thin%struct(n)%map(1,m)
               j=this%ccl_thin%struct(n)%map(2,m)
               k=this%ccl_thin%struct(n)%map(3,m)
               SR_tmp =SR(1,i,j,k)*tangent(1)**2+SR(2,i,j,k)*tangent(2)**2+SR(3,i,j,k)*tangent(3)**2 + &
                  & 2.0_WP*(SR(4,i,j,k)*tangent(1)*tangent(2)+SR(5,i,j,k)*tangent(2)*tangent(3)+SR(6,i,j,k)*tangent(1)*tangent(3))
               lSR(n) = max(lSR(n),abs(SR_tmp))
            end do
         end do
         ! Find the maximum tangential strain rate of each ligament
         call MPI_ALLREDUCE(MPI_IN_PLACE,lSR,1*this%ccl_thin%nstruct,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
         ! Zero out monitoring variables
         this%vof_tf_lig=0.0_WP
         this%np_lig=0
      
   
   
   
   
   
   
   
   
   
         ! if (this%num_old_id.gt.0) then
         ! allocate(overlap(this%ccl_thin%nstruct,this%num_old_id));overlap=0.0_WP
         ! do n=1,this%ccl_thin%nstruct
         !    do m=1,this%ccl_thin%struct(n)%n_
         !       ! Get cell indices
         !       i=this%ccl_thin%struct(n)%map(1,m)
         !       j=this%ccl_thin%struct(n)%map(2,m)
         !       k=this%ccl_thin%struct(n)%map(3,m)
               
         !       if (this%old_id(i,j,k).gt.0) overlap(this%ccl_thin%id(i,j,k),this%old_id(i,j,k)) = overlap(this%ccl_thin%id(i,j,k),this%old_id(i,j,k)) + 1
         !    end do
         ! end do
         ! call MPI_ALLREDUCE(MPI_IN_PLACE, overlap, this%ccl_thin%nstruct*this%num_old_id, MPI_INTEGER, MPI_SUM, this%cfg%comm, ierr)
         ! do n=1, this%ccl_thin%nstruct
         !    if (maxval(overlap(n,:)).gt.0.0_WP.and.this%lig_timers(1,1).gt.VFlo) then
         !       l = maxloc(overlap(n,:),dim=1)
         !       lig_timers_new(n,1) = this%lig_timers(l,1) + dt
         !       lig_timers_new(n,2) = this%lig_timers(l,2)
         !    end if
         ! end do
         ! end if
         ! this%old_id = 0
         ! do n=1,this%ccl_thin%nstruct
         !    do m=1,this%ccl_thin%struct(n)%n_
         !       ! Get cell indices
         !       i=this%ccl_thin%struct(n)%map(1,m)
         !       j=this%ccl_thin%struct(n)%map(2,m)
         !       k=this%ccl_thin%struct(n)%map(3,m)
               
         !       this%old_id(i,j,k)=this%ccl_thin%id(i,j,k)
         !    end do
         ! end do
   
         ! Perform transfer
         do n=1,this%ccl_thin%nstruct
            frem_active = .false.
            if (this%use_film_transfer .and. (lnum(n).ge.9) .and. (lthc(n).le.this%fmin) .and. (lvol(n).gt.0.25*this%cfg%min_meshsize**3) .and. (lper(n).le.0.3)) then
               ! Too close to the end of domain
               !else if (frem(n).gt.0.0_WP) then
               !frem_active = .true.
               ! output to confirm
               lig_film(n) = 1
               if (this%cfg%amRoot) print *, "This is a thin film with min_thickness", lthc(n), "and this is id:", n ,"vol is:", lvol(n) ,"min vol is:", 0.25*this%cfg%min_meshsize**3 ,"num cells is:", lnum(n) ,"lper is:", lper(n)
               ! Assume fd0 across the processor based on the total volume of the film
               if (.not.frem_active) then
                  fd0 =(6.0_WP*Pi*lvol(n)/this%fbvol2dvol)**(1.0_WP/3.0_WP)
               else
                  fd0 =0.0025_WP
               end if
               ! sort cell index based on local film thickness
               if (this%ccl_thin%struct(n)%n_.ge.1) then
                  allocate(sort_id(1:this%ccl_thin%struct(n)%n_)) 
                  allocate(sort_ke(1:this%ccl_thin%struct(n)%n_))
                  ncell_ = this%ccl_thin%struct(n)%n_
                  do m=1,ncell_
                     i=this%ccl_thin%struct(n)%map(1,m)
                     j=this%ccl_thin%struct(n)%map(2,m)
                     k=this%ccl_thin%struct(n)%map(3,m)
                     sort_id(m)=m 
                     sort_ke(m)=this%vf%thickness(i,j,k)
                  end do
                  ! sort based on thickness
                  do ii = 1, ncell_-1
                     do jj = 1, ncell_-ii
                        if (sort_ke(jj).gt.sort_ke(jj+1)) then
                           ! Swap the values
                           tmp_ke = sort_ke(jj)
                           sort_ke(jj) = sort_ke(jj+1)
                           sort_ke(jj+1) = tmp_ke
                           ! Swap the corresponding IDs
                           tmp_id = sort_id(jj)
                           sort_id(jj) = sort_id(jj+1)
                           sort_id(jj+1) = tmp_id
                        end if
                     end do
                  end do
                  Vt=0.0_WP; Vl=0.0_WP
                  np_old=this%lp_spray%np_
                  do m=1,ncell_
                     i=this%ccl_thin%struct(n)%map(1,sort_id(m))
                     j=this%ccl_thin%struct(n)%map(2,sort_id(m))
                     k=this%ccl_thin%struct(n)%map(3,sort_id(m))
                     ! Accumulate 
                     Vl=Vl+this%vf%VF(i,j,k)*this%cfg%vol(i,j,k)
                     if (.not.sampled) then
                        ! Get droplet information based on localized curvature
                        curv_sum=0.0_WP; ncurv=0.0_WP
                        do l=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
                              curv_sum=curv_sum+abs(this%vf%curv2p(l,i,j,k))
                              ncurv=ncurv+1.0_WP
                           end if
                        end do
                        ! call bag_droplet_gamma(this%vf%thickness(i,j,k),2.0_WP*ncurv/curv_sum)
                        ! call bag_droplet_gamma(this%fmin,2.0_WP*ncurv/curv_sum)
                        if (curv_sum.gt.1.0e-14_WP) then
                           call bag_droplet_gamma(this%fmin,ncurv/curv_sum)
                        else
                           call bag_droplet_gamma(this%fmin,0.0_WP)
                        end if
                        Vd = pi/6.0_WP*(min(random_gamma(alpha)*beta*fd0,2.0_WP*frp))**3
                        sampled = .true.
                     end if
                     if (Vl.gt.Vd) then
                        nref=calculateNormal(this%vf%interface_polygon(1,i,j,k))
                        select case (maxloc(abs(nref),1))
                        case (1)
                           tref=normalize([+nref(2),-nref(1),0.0_WP])
                        case (2)
                           tref=normalize([0.0_WP,+nref(3),-nref(2)])
                        case (3)
                           tref=normalize([-nref(3),0.0_WP,+nref(1)])
                        end select
                        sref=cross_product(nref,tref)
   
                        this%lp_spray%np_=this%lp_spray%np_+1
                        ! Make room for new drop
                        call this%lp_spray%resize(this%lp_spray%np_)
                        ! Add the drop
                        if (frem_active) then
                           !this%lp_spray%p(this%lp_spray%np_)%id  =int(7,8)
                        else                                   
                           !this%lp_spray%p(this%lp_spray%np_)%id  =int(6,8)
                           this%lp_spray%p(this%lp_spray%np_)%id  =int(3,8)
                        end if
                        this%lp_spray%p(this%lp_spray%np_)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            
                        this%lp_spray%p(this%lp_spray%np_)%pos =this%vf%Lbary(:,i,j,k)+random_uniform(-0.5_WP*this%cfg%meshsize(i,j,k),0.5_WP*this%cfg%meshsize(i,j,k))*tref+random_uniform(-0.5_WP*this%cfg%meshsize(i,j,k),0.5_WP*this%cfg%meshsize(i,j,k))*sref
                        this%lp_spray%p(this%lp_spray%np_)%vel =this%cfg%get_velocity(pos=this%lp_spray%p(this%lp_spray%np_)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W)    !< Interpolate local cell velocity as drop velocity
                        this%lp_spray%p(this%lp_spray%np_)%ind =this%lp_spray%cfg%get_ijk_global(this%lp_spray%p(this%lp_spray%np_)%pos,[this%lp_spray%cfg%imin,this%lp_spray%cfg%jmin,this%lp_spray%cfg%kmin])    !< Place the drop in the proper cell for the this%lp_spray%cfg
                        this%lp_spray%p(this%lp_spray%np_)%flag=0                                          
                        this%lp_spray%p(this%lp_spray%np_)%dt  =0.0_WP                                     
                        this%lp_spray%p(this%lp_spray%np_)%Acol=0.0_WP                                     
                        this%lp_spray%p(this%lp_spray%np_)%Tcol=0.0_WP  
                        this%lp_spray%p(this%lp_spray%np_)%t   =0.0_WP  
   
                        ! Update tracked volumes
                        Vl=Vl-Vd
                        Vt=Vt+Vd
                        sampled = .false.
   
                        ! Increment monitoring variables
                        this%vof_tf_film=this%vof_tf_film+Vd
                        this%np_film=this%np_film+1
   
                        this%lp_spray%np_new=this%lp_spray%np_new+1
                        this%lp_spray%vp_new=this%lp_spray%vp_new+Vd
                     end if
                     ! Remove liquid in that cell
                     this%vf%VF(i,j,k)=0.0_WP
                  end do
                  deallocate(sort_id,sort_ke)
                  ! If for some reason a film with 0 liquid volume has been tagged, skip it
                  if (Vt.eq.0.0_WP .and. Vl.eq.0.0_WP) cycle
                  ! Based on how many particles were created, decide what to do with left-over volume
                  if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
                     ! Make room for new drop            
                     this%lp_spray%np_=this%lp_spray%np_+1
                     ! Make room for new drop
                     call this%lp_spray%resize(this%lp_spray%np_)
                     ! Add the drop
                     if (frem_active) then
                        !this%lp_spray%p(this%lp_spray%np_)%id  =int(4,8)
                     else                                   
                        !this%lp_spray%p(this%lp_spray%np_)%id  =int(3,8)
                        this%lp_spray%p(this%lp_spray%np_)%id  =int(4,8)
                     end if
                     this%lp_spray%p(this%lp_spray%np_)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            
                     this%lp_spray%p(this%lp_spray%np_)%pos =this%vf%Lbary(:,i,j,k)      
                     this%lp_spray%p(this%lp_spray%np_)%vel =this%cfg%get_velocity(pos=this%lp_spray%p(this%lp_spray%np_)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) !< Interpolate local cell velocity as drop velocity
                     this%lp_spray%p(this%lp_spray%np_)%ind =this%lp_spray%cfg%get_ijk_global(this%lp_spray%p(this%lp_spray%np_)%pos,[this%lp_spray%cfg%imin,this%lp_spray%cfg%jmin,this%lp_spray%cfg%kmin])    !< Place the drop in the proper cell for the this%lp_spray%cfg
                     this%lp_spray%p(this%lp_spray%np_)%flag=0                                          
                     this%lp_spray%p(this%lp_spray%np_)%dt  =0.0_WP                                     
                     this%lp_spray%p(this%lp_spray%np_)%Acol=0.0_WP                                     
                     this%lp_spray%p(this%lp_spray%np_)%Tcol=0.0_WP  
                     this%lp_spray%p(this%lp_spray%np_)%t   =0.0_WP  
   
   
                     ! Increment monitoring variables
                     this%np_film=this%np_film+1
   
                     this%lp_spray%np_new=this%lp_spray%np_new+1
                  else ! Some particles were created, make them all larger
                     do ip=np_old+1,this%lp_spray%np_
                        this%lp_spray%p(ip)%d=this%lp_spray%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
                     end do
                  end if
                  ! Increment monitoring variables
                  this%vof_tf_film=this%vof_tf_film+Vl
                  
                  this%lp_spray%vp_new=this%lp_spray%vp_new+Vl
               end if
            end if
            if(this%use_lig_transfer) then
               at_edge = .false.
               at_edge_total = .false.
               call fit_spline(n,s_info)
               
               if (s_info%flag.lt.1.0_WP) then
                  Lrim = s_info%length
               else
                  Lrim=0.0_WP!llen(n)
               end if
               Vrim=lvol(n)
               if (Lrim .le. VFlo .or. Vrim .le. VFlo) cycle
               minor_radius = 0.0_WP
               l = 0.0_WP
               do m=1,this%ccl_thin%struct(n)%n_
                  ! Get cell indices
                  i=this%ccl_thin%struct(n)%map(1,m)
                  j=this%ccl_thin%struct(n)%map(2,m)
                  k=this%ccl_thin%struct(n)%map(3,m)
                  if (this%recon_type(i,j,k).eq.1) then
                     a_aligned_Cylinder = getAlignedCylinder(this%vf%liquid_gas_interface(i,j,k))
                     minor_radius = minor_radius + sqrt(a_aligned_Cylinder(1))
                     l=l+1.0_WP
                  end if
                  if (i.le.this%cfg%imin+10.or.i.ge.this%cfg%imax-10.or.j.le.this%cfg%jmin+10.or.j.ge.this%cfg%jmax-10.or.k.le.this%cfg%kmin+10.or.k.ge.this%cfg%kmax-10) then
                     at_edge = .true.
                  end if
               end do
               call MPI_ALLREDUCE(at_edge, at_edge_total, 1, MPI_LOGICAL, MPI_LOR, this%cfg%comm, ierr)
               call MPI_ALLREDUCE(MPI_IN_PLACE,minor_radius,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
               call MPI_ALLREDUCE(MPI_IN_PLACE,l,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
               if (minor_radius.ge.VFlo) then
                  if (l .gt. 0) then
                     minor_radius=minor_radius/l
                  else
                     minor_radius=sqrt(Vrim/pi/Lrim)
                  end if
               else
                  minor_radius=sqrt(Vrim/pi/Lrim)
               end if
               if (minor_radius .le. 0.8*sqrt(Vrim/pi/Lrim)) minor_radius = sqrt(Vrim/pi/Lrim)
               if (minor_radius .le. VFlo) cycle
               Trp = 0.0_WP
               ! Drop size method from Kim & Moin (2020)
               nmain=floor(this%dw*Lrim/(twoPi*minor_radius))
   
               ! Calculate breakup time scale based on inviscid RP instability analysis
               ! if (lig_timers_new(n,2).eq.0) then
               !    Trp=2.91258_WP*sqrt(this%fs%rho_l*minor_radius**3/this%fs%sigma)
               !    lig_timers_new(n,2) = Trp
               !    lig_timers_new(n,1) = dt
               ! else
               !    Trp = lig_timers_new(n,2)
               ! end if
               ! ! Calcuate time scale based on maximum local strainrate
               ! if (lSR(n) .le. VFlo) then
               !    Tsr = HUGE(1.0_WP)
               ! else
               !    Tsr=1.0_WP/lSR(n)
               ! end if
   
               ! Only breakup if minimum thickness is reached and there is at least one main drop
               !if (cfg%amRoot) print *, lig_timers_new(n,1), " ", Trp, " ", dt
               call random_number(prob)
               prob = prob * lthc(n)/this%cfg%min_meshsize

               nu = 1.0_WP
               !d_0 = (6.0_WP * lvol(n) / pi)**(1.0_WP/3.0_WP)
               d_0 = lvol(n)**(1.0_WP/3.0_WP)
               d_avg = 0.4_WP * d_0
               d_break = d_avg / 10.0_WP**(nu/6.5)
               !if (this%cfg%amRoot) print *, "n ", n, " d_0 ", d_0, " d_avg ", d_avg, " minor radius ", minor_radius, " break radius ", d_break/2.0_WP
               
               !if (at_edge_total.or.(.not.((prob.le.0.01_WP).and.(Lrim.gt.4*minor_radius).and.(Vrim.ge.0.01*this%cfg%min_meshsize**3).and.(lnum(n).ge.10).and.(minor_radius.ge.0.01*this%cfg%min_meshsize).and.(nmain.ge.1).and.(lthc(n).le.2*this%cfg%min_meshsize).and.(lper(n).gt.0.9)))) cycle!.and.(Trp.le.dt)!.and.(Trp.le.Tsr).and.(lig_timers_new(n,1).ge.Trp)
               if (at_edge_total.or.(.not.((Lrim.gt.4*minor_radius).and.(Vrim.ge.0.01*this%cfg%min_meshsize**3).and.(lnum(n).ge.10).and.(minor_radius.ge.0.01*this%cfg%min_meshsize).and.(nmain.ge.1).and.(minor_radius.le.d_break/2.0_WP).and.(lper(n).gt.0.9)))) cycle
               lig_film(n) = 2
               if (this%cfg%amRoot) then
                  if (this%cfg%amRoot) print *, "This is the min_thickness", lthc(n), "radius:", minor_radius, "lig percentage:", lper(n),"max length:",Lrim,&
                     & "how many cells",lnum(n), "vol:",Lvol(n),"nmain", nmain, "and id:", n !"Trp:", Trp, "Tsr:", Tsr, "Trp/Tsr", Trp/Tsr,
                  
                  nsat=nmain+1
                  diam=(6.0_WP*Vrim/pi/(real(nmain,WP)+this%size_ratio**3*real(nsat,WP)))**(1.0_WP/3.0_WP)

                  filename='spray-all/droplets'
                  open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
                  !allocate(points(3,nsat+nmain))
                  !if (s_info%flag.eq.0) call distribute_on_spline(nsat+nmain,s_info,points)

                  if(allocated(lig_sizes)) deallocate(lig_sizes)
                  call lig_break_gamma(lvol(n),nu,d_avg,lig_sizes)
                  nmain = size(lig_sizes)

                  allocate(points(3,nmain))
                  if (s_info%flag.eq.0) call distribute_on_spline(nmain,s_info,points)

                  !do l=1,nsat+nmain
                  do l=1,nmain
                  ! Only the main processor is in charge of creating droplets
   
                     if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
   
                     this%lp_spray%np_ = this%lp_spray%np_ + 1
                     ! Make room for new drop
                     call this%lp_spray%resize(this%lp_spray%np_)
                     ! Add the drop
                     this%lp_spray%p(this%lp_spray%np_)%id = int(2, 8)
   
                     ! Set the diameter correctly
                     !if (mod(l, 2) .eq. 1) then
                     !   this%lp_spray%p(this%lp_spray%np_)%d = diam * this%size_ratio
                     !else
                     !   this%lp_spray%p(this%lp_spray%np_)%d = diam
                     !end if
                     this%lp_spray%p(this%lp_spray%np_)%d = lig_sizes(l)
         
                     if (s_info%flag.ne.0) then
                        this%lp_spray%p(this%lp_spray%np_)%pos = lpos(n,:)
                     else
                        this%lp_spray%p(this%lp_spray%np_)%pos = points(:,l)
                     end if
                     this%lp_spray%p(this%lp_spray%np_)%vel = lvel(n,:)
         
                     ! Set the remaining particle properties
                     this%lp_spray%p(this%lp_spray%np_)%ind = this%lp_spray%cfg%get_ijk_global(this%lp_spray%p(this%lp_spray%np_)%pos, [this%lp_spray%cfg%imin, this%lp_spray%cfg%jmin, this%lp_spray%cfg%kmin])
                     this%lp_spray%p(this%lp_spray%np_)%flag = 0
                     this%lp_spray%p(this%lp_spray%np_)%dt = 0.0_WP
                     this%lp_spray%p(this%lp_spray%np_)%Acol = 0.0_WP
                     this%lp_spray%p(this%lp_spray%np_)%Tcol = 0.0_WP
                     this%lp_spray%p(this%lp_spray%np_)%t   =0.0_WP  
         
                     ! Output diameter, velocity, and position
                     write(iunit, *) this%lp_spray%p(this%lp_spray%np_)%d, this%lp_spray%p(this%lp_spray%np_)%vel(1), &
                        & this%lp_spray%p(this%lp_spray%np_)%vel(2), this%lp_spray%p(this%lp_spray%np_)%vel(3), &
                        & norm2(this%lp_spray%p(this%lp_spray%np_)%vel), this%lp_spray%p(this%lp_spray%np_)%pos(1), &
                        & this%lp_spray%p(this%lp_spray%np_)%pos(2), this%lp_spray%p(this%lp_spray%np_)%pos(3), this%lp_spray%p(this%lp_spray%np_)%id, 0
                     !print *, "I wrote one particle out of", nsat + nmain
                        print *, "I wrote one particle out of", nmain
                  end do
                  if(allocated(lig_sizes)) deallocate(lig_sizes)
                  ! Close the file
                  close(iunit)
                  ! Increment monitoring variables
                  !this%lp_spray%np_new=this%lp_spray%np_new+nmain+nsat
                  this%lp_spray%np_new=this%lp_spray%np_new+nmain
                  this%lp_spray%vp_new=this%lp_spray%vp_new+lvol(n)
                  !this%np_lig=this%np_lig+nmain+nsat
                  this%np_lig=this%np_lig+nmain
                  this%vof_tf_lig=this%vof_tf_lig+lvol(n)
                  deallocate(points)
               end if
               ! empty out the VF
               do m=1,this%ccl_thin%struct(n)%n_
                  i=this%ccl_thin%struct(n)%map(1,m); j=this%ccl_thin%struct(n)%map(2,m); k=this%ccl_thin%struct(n)%map(3,m)
                  this%vf%VF(i,j,k)=0.0_WP
                  !this%old_id(i,j,k) = 0
               end do
               !lig_timers_new(n,:) = 0.0_WP
            end if
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,lig_film,this%ccl_thin%nstruct,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
   
         ! Gather the number of newly generated particles from each processor due to film burst
         totalnewp = 0
         allocate(plist(0:this%cfg%nproc-1))
         ! Get number of particle generated for each processor
         call MPI_AllGATHER(this%np_film,1,MPI_INTEGER,plist,1,MPI_INTEGER,this%cfg%comm,ierr)
         totalnewp= sum(plist)
         ! If there is any particle generated
         if (totalnewp .gt. 0) then
            ! Transposed allocation for safe, contiguous MPI memory blocks
            allocate(pinfo_(1:max(1,this%np_film), 1:10))
            allocate(pinfo(1:max(1,totalnewp), 1:10))
            allocate(dispels(0:this%cfg%nproc-1))
            
            ! Get info
            count = 0
            do ip = np_start+1, this%lp_spray%np_
               if (this%lp_spray%p(ip)%id .ne. 2) then
                  count = count + 1
                  if (count .le. this%np_film) then
                     pinfo_(count,1) = this%lp_spray%p(ip)%d
                     pinfo_(count,2) = this%lp_spray%p(ip)%vel(1)
                     pinfo_(count,3) = this%lp_spray%p(ip)%vel(2)
                     pinfo_(count,4) = this%lp_spray%p(ip)%vel(3)
                     pinfo_(count,5) = norm2(this%lp_spray%p(ip)%vel)
                     pinfo_(count,6) = this%lp_spray%p(ip)%pos(1)
                     pinfo_(count,7) = this%lp_spray%p(ip)%pos(2)
                     pinfo_(count,8) = this%lp_spray%p(ip)%pos(3)
                     pinfo_(count,9) = real(this%lp_spray%p(ip)%id, WP)
                     pinfo_(count,10) = 0.0_WP
                  end if
               end if
            end do
            
            ! Calculate dispels
            count = 0
            do rank=0,this%cfg%nproc-1
               dispels(rank) = count
               count = count + plist(rank)
            end do
            
            ! Communicate to root
            do i = 1,10
               call MPI_GATHERV(pinfo_(:,i),this%np_film,MPI_REAL_WP,pinfo(:,i),plist,dispels,MPI_REAL_WP,0,this%cfg%comm,ierr)
            end do
            
            !!! Write to droplet list !!!
            if (this%cfg%amRoot)  then
               filename='spray-all/droplets'
               open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
               if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
               do i = 1,totalnewp
                  write(iunit,'(8(f24.16,1x),2(I2,1x))') pinfo(i,1), pinfo(i,2), pinfo(i,3), &
                  & pinfo(i,4), pinfo(i,5), pinfo(i,6), pinfo(i,7), pinfo(i,8), INT(pinfo(i,9)), INT(pinfo(i,10))
               end do
               close(iunit)
            end if
            deallocate(pinfo_, pinfo, dispels)
            ! Synchronize VF fields
            call this%cfg%sync(this%vf%VF)
            call this%vf%clean_irl_and_band()
            ! Synchronize particles
            call this%lp_spray%sync()
            ! Integrate monitoring variables 
            call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_tf_film,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_film    ,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
         end if
   
         ! Synchronize VF fields
         call this%cfg%sync(this%vf%VF)
         call this%vf%clean_irl_and_band()
         ! Synchronize particles
         call this%lp_spray%sync()
         ! Integrate monitoring variables
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_tf_lig,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_lig    ,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)

         ! survivor_count = 0
         ! do n = 1, this%ccl_thin%nstruct
         !    if (lig_film(n).eq.2 .and. lig_timers_new(n, 1) .gt. VFlo) then
         !       survivor_count = survivor_count + 1
         !    end if
         ! end do

         ! if (survivor_count.gt.0 .and. survivor_count.lt.this%ccl_thin%nstruct) then
         !    allocate(packed_lig_timers(survivor_count, 2))
         !    allocate(survivor_map(this%ccl_thin%nstruct))
         !    survivor_map = 0
            
         !    current_id = 0
         !    do n = 1, this%ccl_thin%nstruct
         !       if (lig_film(n).eq.2 .and. lig_timers_new(n, 1) .gt. VFlo) then
         !          current_id = current_id + 1
         !          packed_lig_timers(current_id, :) = lig_timers_new(n, :)
         !          survivor_map(n) = current_id
         !       end if
         !    end do
            
         !    do k=this%cfg%kmin_,this%cfg%kmax_
         !       do j=this%cfg%jmin_,this%cfg%jmax_
         !          do i=this%cfg%imin_,this%cfg%imax_
         !             if (this%old_id(i,j,k) .gt. 0) then
         !                this%old_id(i,j,k) = survivor_map(this%old_id(i,j,k))
         !             end if
         !          end do
         !       end do
         !    end do
         !    deallocate(survivor_map)
            
         !    if (allocated(this%lig_timers)) deallocate(this%lig_timers)
         !    allocate(this%lig_timers(survivor_count, 2))
         !    this%lig_timers = packed_lig_timers
         !    this%num_old_id = survivor_count
         !    deallocate(packed_lig_timers)
         ! else if (survivor_count.gt.0) then
         !    if (allocated(this%lig_timers)) deallocate(this%lig_timers)
         !    allocate(this%lig_timers(this%ccl_thin%nstruct,2))
         !    this%lig_timers = lig_timers_new
         !    this%num_old_id = this%ccl_thin%nstruct
         ! else
         !    this%num_old_id = 0
         !    if (allocated(this%lig_timers)) deallocate(this%lig_timers)
         !    allocate(this%lig_timers(1,1)); this%lig_timers = 0.0_WP
         ! end if
      else
         ! this%num_old_id = 0
         ! if (allocated(this%lig_timers)) deallocate(this%lig_timers)
         ! allocate(this%lig_timers(1,1)); this%lig_timers = 0.0_WP
      end if
      if (allocated(lig_film)) deallocate(lig_film)
      if (allocated(lvol)) deallocate(lvol)
      if (allocated(lthc)) deallocate(lthc)
      if (allocated(llen)) deallocate(llen)
      if (allocated(lnum)) deallocate(lnum)
      if (allocated(lper)) deallocate(lper)
      if (allocated(lpos)) deallocate(lpos)
      if (allocated(lvel)) deallocate(lvel)
      if (allocated(lmoi)) deallocate(lmoi)
      if (allocated(lSR)) deallocate(lSR)
      if (allocated(xmin)) deallocate(xmin, xmax)
      if (allocated(ymin)) deallocate(ymin, ymax)
      if (allocated(zmin)) deallocate(zmin, zmax)
      if (allocated(SR)) deallocate(SR)
      if (allocated(lig_timers_new)) deallocate(lig_timers_new)
      if (allocated(local_thickness)) deallocate(local_thickness)
      if (allocated(local_struct_type)) deallocate(local_struct_type)
      if (allocated(overlap)) deallocate(overlap)
   
   contains
      ! Calculate thickness and struct_type based on moment of inertia
      subroutine get_structinfo()
         implicit none
         real(WP) :: tmpvol,tmparea
         real(WP), dimension(1:3) :: tmpxvol, tmpL
         integer :: nneigh_moi, nneigh_thickness
         real(WP) :: x1,x2,phi
         nneigh_moi=2; nneigh_thickness=3
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (this%vf%VF(i,j,k) .le. VFlo) cycle
                  ! calculate thickness
                  tmpvol=0.0_WP; tmparea=0.0_WP
                  do kk = k-nneigh_thickness,k+nneigh_thickness
                     do jj = j-nneigh_thickness,j+nneigh_thickness
                        do ii = i-nneigh_thickness,i+nneigh_thickness
                           tmpvol = tmpvol + this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
                           if (this%recon_type(ii,jj,kk).eq.2.or.this%recon_type(ii,jj,kk).eq.3.or.this%recon_type(ii,jj,kk).eq.0) then
                              tmparea = tmparea + this%vf%SD(ii,jj,kk)*this%cfg%vol(i,j,k)
                           else
                              tmparea = tmparea + this%vf%SD(ii,jj,kk)*this%cfg%vol(i,j,k)*(2.0/sqrt(pi))
                           end if
                        end do
                     end do
                  end do
                  ! Calculate thickness
                  if (this%vf%VF(i,j,k).le.VFlo) then
                     local_thickness(i,j,k) = 0.0_WP
                  else if (tmparea .gt. 0.0_WP) then
                     local_thickness(i,j,k) = 2.0_WP*tmpvol/(tmparea+tiny(1.0_WP))
                  else
                     local_thickness(i,j,k) = 3.5_WP*this%cfg%min_meshsize
                  end if
   
                  if (local_thickness(i,j,k).lt.1.0_WP*this%cfg%min_meshsize) then
                     ! Calculate moi
                     tmpvol=0.0_WP; tmpxvol=0.0_WP; A=0.0_WP
                     ! First pass to accumulate volume, surface area, and position
                     do kk = k-nneigh_moi,k+nneigh_moi
                        do jj = j-nneigh_moi,j+nneigh_moi
                           do ii = i-nneigh_moi,i+nneigh_moi
                              tmpvol = tmpvol + this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
                              tmpxvol = tmpxvol + this%vf%Lbary(:,ii,jj,kk)*this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
                           end do
                        end do
                     end do
                     ! Second pass to accumulate moment of inertia
                     if (tmpvol.gt.1.0e-14_WP) tmpxvol = tmpxvol/tmpvol
                     do kk = k-nneigh_moi,k+nneigh_moi
                        do jj = j-nneigh_moi,j+nneigh_moi
                           do ii = i-nneigh_moi,i+nneigh_moi
                              ! Location of film node
                              tmpL = this%vf%Lbary(:,ii,jj,kk) - tmpxvol
                              A(1,1)=A(1,1)+this%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*(tmpL(2)**2+tmpL(3)**2)
                              A(2,2)=A(2,2)+this%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*(tmpL(1)**2+tmpL(3)**2)
                              A(3,3)=A(3,3)+this%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*(tmpL(1)**2+tmpL(2)**2)
                              A(1,2)=A(1,2)-this%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*tmpL(1)*tmpL(2)
                              A(1,3)=A(1,3)-this%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*tmpL(1)*tmpL(3)
                              A(2,3)=A(2,3)-this%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*tmpL(2)*tmpL(3)
                           end do
                        end do
                     end do
                     A(2,1) = A(1,2)
                     A(3,1) = A(1,3)
                     A(3,2) = A(2,3)
         
                     x1 = A(1,1)**2+A(2,2)**2+A(3,3)**2-A(1,1)*A(2,2)-A(1,1)*A(3,3)-A(2,2)*A(3,3)+3*(A(1,2)**2+A(1,3)**2+A(2,3)**2)
                     x2 = -(2*A(1,1)-A(2,2)-A(3,3))*(2*A(2,2)-A(1,1)-A(3,3))*(2*A(3,3)-A(1,1)-A(2,2))+9.0_WP*((2*A(3,3)-A(1,1)-A(2,2))*A(1,2)**2+(2*A(2,2)-A(1,1)-A(3,3))*A(1,3)**2+(2*A(1,1)-A(2,2)-A(3,3))*A(2,3)**2)-54.0_WP*A(1,2)*A(1,3)*A(2,3)
         
                     phi = atan2(sqrt(max(0.0_WP, 4*x1**3 - x2**2)), x2)
         
                     d(1) = (A(1,1)+A(2,2)+A(3,3)-2*sqrt(max(0.0_WP, x1))*cos(phi/3.0_WP))/3.0_WP
                     d(2) = (A(1,1)+A(2,2)+A(3,3)+2*sqrt(max(0.0_WP, x1))*cos((phi+pi)/3.0_WP))/3.0_WP
                     d(3) = (A(1,1)+A(2,2)+A(3,3)+2*sqrt(max(0.0_WP, x1))*cos((phi-pi)/3.0_WP))/3.0_WP  
                     ! Calculate local struct type
                     d=max(0.0_WP,d)
                     if ((d(3).lt.1.5_WP*d(2)).and.(d(2).gt.1.5_WP*d(1)) ) then
                        local_struct_type(i,j,k) = 1
                     else if ( (d(3) .gt. 1.5_WP * d(2)) .and. (d(2) .lt. 1.5_WP * d(1)) ) then
                        local_struct_type(i,j,k) = 2
                     end if
                  end if
               end do
            end do
         end do
         call this%cfg%sync(local_thickness)
         call this%cfg%sync(local_struct_type)
      end subroutine get_structinfo
   
      subroutine bag_droplet_gamma(h,R)
         implicit none
         real(WP), intent(in) :: h,R
         real(WP) :: Utc,ac,b,dr,ds,Oh
         real(WP) :: mean, stdev
         ! assert h,R != 0
         ! Retraction speed
         Utc=0.0_WP
         if (h.gt.1.0e-14_WP) Utc=sqrt(2.0_WP*this%fs%sigma/this%fs%rho_l/h)
         ! Centripetal acceleration
         ac=0.0_WP
         if (R.gt.1.0e-14_WP) ac=Utc**2/R
         ! Rim diameter
         b=0.0_WP
         if (ac.gt.1.0e-14_WP) b=sqrt(this%fs%sigma/this%fs%rho_l/ac)
         ! RP droplet diameter
         frp=1.89_WP*b
         ! Rim Ohnesorge number
         Oh=0.0_WP
         if (b.gt.1.0e-14_WP) Oh=this%fs%visc_l/sqrt(this%fs%rho_l*b*this%fs%sigma)
         ! Satellite droplet diameter
         ds=frp/sqrt(2.0_WP+3.0_WP*Oh/sqrt(2.0_WP))
         ! Mean and standard deviation of diameter of all modes, normalized by drop diameter
         mean=0.25_WP*(h+b+frp+ds)/fd0
         stdev=sqrt(0.25_WP*sum(([h,b,frp,ds]/fd0-mean)**2))
         ! Gamma distribution parameters
         alpha = 0
         beta = 0
         if(stdev.gt.1.0e-14_WP) alpha=(mean/stdev)**2
         if(mean.gt.1.0e-14_WP) beta=stdev**2/mean
      end subroutine bag_droplet_gamma

      subroutine lig_break_gamma(vol,n,d_mean,sizes)
         use mathtools, only: pi
         use random, only: random_gamma
         use vfs_data_class, only: VFlo,VFhi
   
         real(WP), intent(in) :: vol, d_mean, n
         real(WP) :: new_vol, vol_ratio
         real(WP), dimension(:), allocatable, intent(out) :: sizes(:)
         integer :: i, num_drops
   
         new_vol = 0.0_WP
         vol_ratio = 0.0_WP
         num_drops = ceiling(vol / (4.0_WP/3.0_WP * pi * (d_mean/2.0_WP)**3))
         !print *, "vol ", vol, " d_mean ", d_mean, " num_drops ", num_drops
         if(allocated(sizes)) deallocate(sizes)
         allocate(sizes(num_drops))
   
         do i=1,num_drops
            sizes(i) = d_mean*random_gamma(n)/n
            new_vol = new_vol + 4.0_WP/3.0_WP * pi * (sizes(i)/2.0_WP)**3
         end do
   
         if (new_vol.ge.VFlo) then
            vol_ratio = (vol/new_vol)**(1.0_WP/3.0_WP)
         else
            vol_ratio = 1.0_WP
         end if
   
         do i=1,num_drops
            sizes(i) = sizes(i) * vol_ratio
         end do
      end subroutine lig_break_gamma
   
      !> Function that identifies cells that need a label
      logical function make_label(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         if ((local_struct_type(i,j,k).eq.1).and.(this%vf%VF(i,j,k).gt.VFlo).and.(local_thickness(i,j,k).lt.2.0*this%cfg%min_meshsize)) then
            make_label=.true.
         else if ((local_struct_type(i,j,k).ne.1).and.(this%vf%VF(i,j,k).gt.VFlo).and.(local_thickness(i,j,k).lt.1.0*this%cfg%min_meshsize))then
            make_label=.true.
         else
            make_label=.false.
         end if
      end function make_label
   
      !> Function that identifies if cell pairs have same label
      logical function same_label(i1,j1,k1,i2,j2,k2)
         implicit none
         integer, intent(in) :: i1,j1,k1,i2,j2,k2
         same_label=.true.
      end function same_label
   
      subroutine fit_spline(n,s_info)
         use fitpack_core, only: curfit, splev
         integer, intent(in) :: n
         type(spline_info), intent(out) :: s_info
         integer :: m, local_point_count, num_procs, total_points, i, j, k, ier, nest_max, lwrk, unique_count
         real(WP), dimension(3) :: end_p
         real(WP), dimension(:,:), allocatable :: local_points
         real(WP), dimension(:,:), allocatable :: points, unique_points, sorted_points
         integer, dimension(:), allocatable :: recv_counts
         integer, dimension(:), allocatable :: displacements
         real(WP), dimension(:), allocatable :: t_param
         real(WP) :: s, fp_x, fp_y, fp_z
         real(WP), dimension(:), allocatable :: weights
         logical, dimension(:), allocatable :: is_used
         real(WP) :: min_dist_sq, dist_sq
         real(WP), dimension(3) :: last_sorted_point
         integer :: start_point_idx, best_idx
   
         real(WP), dimension(:), allocatable :: wrk
         integer, dimension(:), allocatable :: iwrk
         real(WP) :: tolerance
   
         integer :: num_eval
         real(WP), dimension(:), allocatable :: t_eval, x_eval, y_eval, z_eval
         real(WP) :: step, t_current
         integer :: e_flag
         integer, allocatable :: cluster_id(:)
         integer, allocatable :: cluster_counts(:)
         character(len=20) :: my_string

         s_info%length = 0.0_WP
         s_info%flag = 0
   
         allocate(local_points(3, this%ccl_thin%struct(n)%n_)); local_points = 0.0_WP
         write (my_string, '(i0)') n
         do m=1,this%ccl_thin%struct(n)%n_
            i=this%ccl_thin%struct(n)%map(1,m)
            j=this%ccl_thin%struct(n)%map(2,m)
            k=this%ccl_thin%struct(n)%map(3,m)
   
            if (this%recon_type(i,j,k).eq.1) then
               local_points(:,m) = getDatum(this%vf%liquid_gas_interface(i,j,k))
               local_points(1,m) = local_points(1,m)-this%ccl_thin%struct(n)%per(1)*this%cfg%xL
               local_points(2,m) = local_points(2,m)-this%ccl_thin%struct(n)%per(2)*this%cfg%yL
               local_points(3,m) = local_points(3,m)-this%ccl_thin%struct(n)%per(3)*this%cfg%zL
            else
               local_points(:,m) = this%vf%LBARY(:,i,j,k)
               local_points(1,m) = local_points(1,m)-this%ccl_thin%struct(n)%per(1)*this%cfg%xL
               local_points(2,m) = local_points(2,m)-this%ccl_thin%struct(n)%per(2)*this%cfg%yL
               local_points(3,m) = local_points(3,m)-this%ccl_thin%struct(n)%per(3)*this%cfg%zL
            end if
         end do
   
         num_procs = this%cfg%nproc
   
         allocate(recv_counts(num_procs))
         recv_counts = 0
   
         call MPI_GATHER(this%ccl_thin%struct(n)%n_, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, 0, this%cfg%comm, ierr)
   
         if (this%cfg%amRoot) then
            total_points = sum(recv_counts)
            allocate(points(3, total_points))
            allocate(is_used(total_points))
            allocate(displacements(num_procs))
            if (total_points.gt.0) then
               displacements(1) = 0
               do i=2,num_procs
                  displacements(i) = displacements(i-1) + recv_counts(i-1)
               end do
               is_used = .false.
            end if
            call MPI_GATHERV(local_points, this%ccl_thin%struct(n)%n_*3, MPI_REAL_WP, points, recv_counts*3, displacements*3, MPI_REAL_WP, 0, this%cfg%comm, ierr)
         else
            call MPI_GATHERV(local_points, this%ccl_thin%struct(n)%n_*3, MPI_REAL_WP, local_points, recv_counts, recv_counts, MPI_REAL_WP, 0, this%cfg%comm, ierr)
         end if
   
         if (this%cfg%amRoot) then 
            tolerance = (4*this%cfg%min_meshsize)**2
            if (total_points .gt. 4) then
               s_info%flag = 0
               end_p = points(:,1)
               start_point_idx = 1
   
               allocate(cluster_id(total_points))
               cluster_id = 0
               unique_count = 0
   
               do m = 1, total_points
                  if (cluster_id(m) .eq. 0) then
                     unique_count = unique_count + 1
                     cluster_id(m) = unique_count
                     do i = m + 1, total_points
                        if (cluster_id(i) .eq. 0) then
                           if (sum((points(:, m) - points(:, i))**2) .le. tolerance) then
                              cluster_id(i) = unique_count
                           end if
                        end if
                     end do
                  end if
               end do
               if (unique_count .gt. 0) then
                  allocate(unique_points(3, unique_count), cluster_counts(unique_count))
                  unique_points = 0.0_WP
                  cluster_counts = 0
                  do m = 1, total_points
                     unique_points(:, cluster_id(m)) = unique_points(:, cluster_id(m)) + points(:, m)
                     cluster_counts(cluster_id(m)) = cluster_counts(cluster_id(m)) + 1
                  end do
                  do m = 1, unique_count
                     if (cluster_counts(m) .gt. 0) then
                        unique_points(:, m) = unique_points(:, m) / real(cluster_counts(m))
                     end if
                  end do
               end if
               deallocate(points)
               allocate(points(3, unique_count))
               allocate(sorted_points(3, total_points))
               points = unique_points
               total_points = unique_count
               deallocate(cluster_id, unique_points, cluster_counts)
   
               if (total_points .gt. 4) then
                  call order_points(total_points, points, unique_count, sorted_points)
                  total_points = unique_count
                  deallocate(points)
                  allocate(points(3, total_points))
                  points = sorted_points
                  deallocate(sorted_points)
                  allocate(t_param(total_points))
                  allocate(weights(total_points))
                  t_param(1) = 0.0_WP
                  do m = 2, total_points
                     t_param(m) = t_param(m-1)+sqrt(sqrt((points(1,m)-points(1,m-1))**2+(points(2,m)-points(2,m-1))**2+(points(3,m)-points(3,m-1))**2))
                  end do
                  if (t_param(total_points) .gt. VFlo) then
                     t_param = t_param/t_param(total_points)
                  else
                     s_info%flag = 1
                  end if
                  s = real(total_points, WP) * (this%cfg%min_meshsize)**2
                  weights = 1.0_WP
                  k = 3
                  nest_max = max(total_points+k+1, 2*k+3)
                  lwrk = total_points * (k + 1) + nest_max * (7 + 3 * k)
      
                  if (allocated(s_info%t_knots_x)) deallocate(s_info%t_knots_x)
                  if (allocated(s_info%c_coeffs_x)) deallocate(s_info%c_coeffs_x)
                  if (allocated(s_info%t_knots_y)) deallocate(s_info%t_knots_y)
                  if (allocated(s_info%c_coeffs_y)) deallocate(s_info%c_coeffs_y)
                  if (allocated(s_info%t_knots_z)) deallocate(s_info%t_knots_z)
                  if (allocated(s_info%c_coeffs_z)) deallocate(s_info%c_coeffs_z)
      
                  allocate(wrk(lwrk))
                  allocate(iwrk(nest_max))
                  allocate(s_info%t_knots_x(nest_max))
                  allocate(s_info%t_knots_y(nest_max))
                  allocate(s_info%t_knots_z(nest_max))
                  allocate(s_info%c_coeffs_x(nest_max))
                  allocate(s_info%c_coeffs_y(nest_max))
                  allocate(s_info%c_coeffs_z(nest_max))
      
                  call curfit(iopt=0, m=total_points, x=t_param, y=points(1,:), w=weights, &
                  xb=t_param(1), xe=t_param(total_points), k=k, s=s, nest=nest_max, &
                  n=s_info%n_knots_x, t=s_info%t_knots_x, c=s_info%c_coeffs_x, fp=fp_x, &
                  wrk=wrk, lwrk=lwrk, iwrk=iwrk, ier=ier)
                  if (ier .gt. 0) s_info%flag = 1!print *, "Error in CURFIT for X: ", ier
      
                  call curfit(iopt=0, m=total_points, x=t_param, y=points(2,:), w=weights, &
                        xb=t_param(1), xe=t_param(total_points), k=k, s=s, nest=nest_max, &
                        n=s_info%n_knots_y, t=s_info%t_knots_y, c=s_info%c_coeffs_y, fp=fp_y, &
                        wrk=wrk, lwrk=lwrk, iwrk=iwrk, ier=ier)
                  if (ier .gt. 0) s_info%flag = 1!print *, "Error in CURFIT for Y: ", ier
      
                  call curfit(iopt=0, m=total_points, x=t_param, y=points(3,:), w=weights, &
                        xb=t_param(1), xe=t_param(total_points), k=k, s=s, nest=nest_max, &
                        n=s_info%n_knots_z, t=s_info%t_knots_z, c=s_info%c_coeffs_z, fp=fp_z, &
                        wrk=wrk, lwrk=lwrk, iwrk=iwrk, ier=ier)
                  if (ier .gt. 0) s_info%flag = 1!print *, "Error in CURFIT for Z: ", ier
      
                  if (s_info%flag.eq.0) then
                     num_eval = total_points*10
                     e_flag = 0
                     allocate(t_eval(num_eval))
                     allocate(x_eval(num_eval))
                     allocate(y_eval(num_eval))
                     allocate(z_eval(num_eval))
                     
                     step = t_param(total_points) / real(num_eval - 1, WP)
                     do m = 1, num_eval
                        t_eval(m) = real(m-1, WP) * step
                     end do
                     
                     call splev(s_info%t_knots_x, s_info%n_knots_x, s_info%c_coeffs_x, k, t_eval, x_eval, num_eval, e_flag, ier)
                     if (ier .ne. 0) print *, "Error in SPLEV for X: ", ier
                     
                     call splev(s_info%t_knots_y, s_info%n_knots_y, s_info%c_coeffs_y, k, t_eval, y_eval, num_eval, e_flag, ier)
                     if (ier .ne. 0) print *, "Error in SPLEV for Y: ", ier
                     
                     call splev(s_info%t_knots_z, s_info%n_knots_z, s_info%c_coeffs_z, k, t_eval, z_eval, num_eval, e_flag, ier)
                     if (ier .ne. 0) print *, "Error in SPLEV for Z: ", ier
      
                     s_info%length = 0.0_WP
                     do i = 2, num_eval
                        s_info%length = s_info%length + sqrt((x_eval(i)-x_eval(i-1))**2+(y_eval(i)-y_eval(i-1))**2+(z_eval(i)-z_eval(i-1))**2)
                     end do
                  end if
               else
                  s_info%flag = 1
               end if
            else
               s_info%flag = 1
            end if
         end if
         call MPI_BCAST(s_info%flag, 1, MPI_INTEGER, 0, this%cfg%comm, ierr)
         call MPI_BCAST(s_info%length, 1, MPI_REAL_WP, 0, this%cfg%comm, ierr)
         if (allocated(local_points)) deallocate(local_points)
         if (allocated(recv_counts)) deallocate(recv_counts)
         if (allocated(points)) deallocate(points)
         if (allocated(is_used)) deallocate(is_used)
         if (allocated(displacements)) deallocate(displacements)
         if (allocated(t_param)) deallocate(t_param)
         if (allocated(weights)) deallocate(weights)
         if (allocated(wrk)) deallocate(wrk)
         if (allocated(iwrk)) deallocate(iwrk)
         if (allocated(t_eval)) deallocate(t_eval)
         if (allocated(x_eval)) deallocate(x_eval)
         if (allocated(y_eval)) deallocate(y_eval)
         if (allocated(z_eval)) deallocate(z_eval)
      end subroutine fit_spline
   
      subroutine distribute_on_spline(n_part,s_info,points)
         use fitpack_core, only: splev
         integer, intent(in) :: n_part
         type(spline_info), intent(in) :: s_info
         real(WP), dimension(:,:), intent(out) :: points
         integer :: k, ier, n_map, idx, i, m
         real(WP), dimension(:), allocatable :: t_map, x_map, y_map, z_map, dist_map
         real(WP) :: step, total_len, target_dist, step_map, frac
         integer :: e_flag

         if (this%cfg%amRoot) then 
            e_flag = 0
            k = 3
            n_map = 2000

            allocate(t_map(n_map), x_map(n_map), y_map(n_map), z_map(n_map), dist_map(n_map))

            step_map = 1.0_WP / real(n_map - 1, WP)
            do i = 1, n_map
               t_map(i) = real(i-1, WP) * step_map
            end do

            call splev(s_info%t_knots_x, s_info%n_knots_x, s_info%c_coeffs_x, k, t_map, x_map, n_map, e_flag, ier)
            if (ier .ne. 0) print *, "Error in SPLEV for X: ", ier
            call splev(s_info%t_knots_y, s_info%n_knots_y, s_info%c_coeffs_y, k, t_map, y_map, n_map, e_flag, ier)
            if (ier .ne. 0) print *, "Error in SPLEV for Y: ", ier
            call splev(s_info%t_knots_z, s_info%n_knots_z, s_info%c_coeffs_z, k, t_map, z_map, n_map, e_flag, ier)
            if (ier .ne. 0) print *, "Error in SPLEV for Z: ", ier

            dist_map(1) = 0.0_WP
            do i = 2, n_map
               dist_map(i) = dist_map(i-1) + sqrt((x_map(i)-x_map(i-1))**2 + (y_map(i)-y_map(i-1))**2 + (z_map(i)-z_map(i-1))**2)
            end do
            total_len = dist_map(n_map)
            step = total_len / real(n_part + 1, WP)
            
            do m = 1, n_part
               target_dist = real(m, WP) * step
               idx = 1
               do i = 1, n_map-1
                  if (dist_map(i+1) .ge. target_dist) then
                     idx = i
                     exit
                  end if
               end do
               
               if (abs(dist_map(idx+1) - dist_map(idx)) .gt. 1.0e-12_WP) then
                   frac = (target_dist - dist_map(idx)) / (dist_map(idx+1) - dist_map(idx))
               else
                   frac = 0.0_WP
               end if
               
               points(1,m) = x_map(idx) + frac * (x_map(idx+1) - x_map(idx))
               points(2,m) = y_map(idx) + frac * (y_map(idx+1) - y_map(idx))
               points(3,m) = z_map(idx) + frac * (z_map(idx+1) - z_map(idx))
   
               if (this%cfg%xper.and.points(1,m).lt.this%cfg%x(this%cfg%imin)) points(1,m)=points(1,m)+this%cfg%xL
               if (this%cfg%yper.and.points(2,m).lt.this%cfg%y(this%cfg%jmin)) points(2,m)=points(2,m)+this%cfg%yL
               if (this%cfg%zper.and.points(3,m).lt.this%cfg%z(this%cfg%kmin)) points(3,m)=points(3,m)+this%cfg%zL
            end do

            if (allocated(t_map)) deallocate(t_map, x_map, y_map, z_map, dist_map)
         end if
      end subroutine distribute_on_spline

      subroutine order_points(total_points, points, total_ordered_points, ordered_points)
         integer, intent(inout) :: total_points
         real(WP), intent(in) :: points(3, total_points)
         integer, intent(out) :: total_ordered_points
         real(WP), allocatable, intent(out) :: ordered_points(:,:)
         
         logical, allocatable :: in_tree(:)
         real(WP), allocatable :: min_tree_dist(:)
         integer, allocatable :: parent(:)
         integer :: i, j, new_node
         real(WP) :: dist, current_dist
         
         integer, allocatable :: degree(:), node_offset(:), adjacency_list(:), local_offset(:)
         
         integer, allocatable :: queue(:), path(:)
         logical, allocatable :: visited(:)
         integer :: q_head, q_tail, node_a, node_b

         if (total_points .le. 1) then
            total_ordered_points = total_points
            allocate(ordered_points(3, max(1, total_points)))
            if (total_points .eq. 1) ordered_points(:,1) = points(:,1)
            return
         end if

         allocate(parent(total_points), in_tree(total_points), min_tree_dist(total_points))
         in_tree = .false.; parent = 0; min_tree_dist = huge(1.0_WP)
         min_tree_dist(1) = 0.0_WP
         
         do i = 1, total_points
            current_dist = huge(1.0_WP); new_node = -1
            do j = 1, total_points
               if (.not. in_tree(j) .and. min_tree_dist(j) .lt. current_dist) then
                  current_dist = min_tree_dist(j); new_node = j
               end if
            end do
            
            if (new_node .eq. -1) exit
            in_tree(new_node) = .true.
            
            do j = 1, total_points
               if (.not. in_tree(j)) then
                  dist = sum((points(:, new_node) - points(:, j))**2)
                  if (dist .lt. min_tree_dist(j)) then
                     min_tree_dist(j) = dist; parent(j) = new_node
                  end if
               end if
            end do
         end do
         deallocate(in_tree, min_tree_dist)

         allocate(degree(total_points), node_offset(total_points+1), local_offset(total_points), adjacency_list(2*total_points - 2))
         degree = 0
         
         do i = 2, total_points
            j = parent(i)
            if (j .gt. 0) then
               degree(i) = degree(i) + 1
               degree(j) = degree(j) + 1
            end if
         end do
         
         node_offset(1) = 1
         do i = 1, total_points
            node_offset(i+1) = node_offset(i) + degree(i)
            local_offset(i) = node_offset(i)
         end do
         
         do i = 2, total_points
            j = parent(i)
            if (j .gt. 0) then
               adjacency_list(local_offset(i)) = j; local_offset(i) = local_offset(i) + 1
               adjacency_list(local_offset(j)) = i; local_offset(j) = local_offset(j) + 1
            end if
         end do
         deallocate(degree, local_offset, parent)

         allocate(queue(total_points), visited(total_points), path(total_points))
         visited = .false.
         
         q_head = 1; q_tail = 1; queue(q_head) = 1; visited(1) = .true.
         node_a = 1
         
         do while (q_head .le. q_tail)
            new_node = queue(q_head); q_head = q_head + 1
            node_a = new_node 
            do i = node_offset(new_node), node_offset(new_node+1) - 1
               j = adjacency_list(i)
               if (.not. visited(j)) then
                  visited(j) = .true.
                  q_tail = q_tail + 1; queue(q_tail) = j
               end if
            end do
         end do

         visited = .false.; path = 0
         q_head = 1; q_tail = 1; queue(q_head) = node_a; visited(node_a) = .true.
         node_b = node_a
         
         do while (q_head .le. q_tail)
            new_node = queue(q_head); q_head = q_head + 1
            node_b = new_node 
            do i = node_offset(new_node), node_offset(new_node+1) - 1
               j = adjacency_list(i)
               if (.not. visited(j)) then
                  visited(j) = .true.; path(j) = new_node
                  q_tail = q_tail + 1; queue(q_tail) = j
               end if
            end do
         end do

         total_ordered_points = 0; new_node = node_b
         do while (new_node .ne. 0)
            total_ordered_points = total_ordered_points + 1
            new_node = path(new_node)
         end do

         allocate(ordered_points(3, total_ordered_points))
         new_node = node_b; i = 1
         do while (new_node .ne. 0)
            ordered_points(:, i) = points(:, new_node)
            i = i + 1
            new_node = path(new_node)
         end do

         deallocate(node_offset, adjacency_list, queue, visited, path)
      end subroutine order_points

   end subroutine transfer_thin_features

   !> Secondary breakup of drops
   subroutine secondary_break(this)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: pi,normalize,cross_product
      use messager,  only: die
      use vfs_data_class, only: VFlo,VFhi
      use random, only: random_normal
      class(detection), intent(inout) :: this

      integer :: m,n,i,j,k,num_part,iunit,ierr
      integer :: indx(3), indx2(3)
      real(WP) :: r_cr,r,We_cr,We,u_rel_sq,t_bu,t,mean,stdev,vol,mean_vol,new_vol,vol_ratio,num1,num2,num3,a1,a2
      real(WP) :: shift(3), vel_shift(3), old_pos(3), old_vel(3)
      real(WP), allocatable :: new_sizes(:)

      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      We_cr = 6.0_WP
      do m=1,this%lp_spray%np_
         t = this%lp_spray%p(m)%t
         r = this%lp_spray%p(m)%d/2.0_WP
         indx = this%cfg%get_ijk_global(this%lp_spray%p(m)%pos,[this%lp_spray%cfg%imin,this%lp_spray%cfg%jmin,this%lp_spray%cfg%kmin])
         i = indx(1); j = indx(2); k = indx(3)
         u_rel_sq = (this%Ui(i,j,k) - this%lp_spray%p(m)%vel(1))**2 + (this%Vi(i,j,k) - this%lp_spray%p(m)%vel(2))**2 + (this%Wi(i,j,k) - this%lp_spray%p(m)%vel(3))**2
         !u_rel_sq = (We_cr*this%fs%sigma) / (this%fs%rho_g*2.0_WP*this%cfg%min_meshsize)
         if (u_rel_sq.ge.VFlo .and. this%fs%rho_g.ge.VFlo) then
            r_cr = (We_cr*this%fs%sigma) / (this%fs%rho_g*u_rel_sq)
         else
            r_cr = HUGE(1.0_WP)
         end if
         if (this%fs%sigma.ge.VFlo) then
            We = (this%fs%rho_g*u_rel_sq*r) / this%fs%sigma
         else
            We = HUGE(1.0_WP)
         end if
         if (u_rel_sq.ge.VFlo .and. this%fs%rho_g.ge.VFlo) then
            t_bu = sqrt(1.0_WP/3.0_WP) * sqrt(this%fs%rho_l/this%fs%rho_g) * r/sqrt(u_rel_sq)
         else
            t_bu = HUGE(1.0_WP)
         end if

         !if (r.ge.this%cfg%min_meshsize) print *,  "r ", r, " r_cr ", r_cr, " We ", We, " We_cr ", We_cr, " t ", t, " t_bu ", t_bu, " u_rel_sq ", u_rel_sq, "  u_rel_sq_forced ", (We_cr*this%fs%sigma) / (this%fs%rho_g*2.0_WP*this%cfg%min_meshsize)
         !t_bu = 0.1_WP
         !if (r.ge.r_cr) print *, "r ", r, " r_cr ", r_cr, " We ", We, " We_cr ", We_cr, " t ", t, " t_bu ", t_bu, " u_rel_sq ", u_rel_sq
         if (r.ge.r_cr .and. We.ge.We_cr .and. t.ge.t_bu) then
            a1 = 0.6_WP*log(We_cr/We)
            a2 = -a1*We
            mean = log(r) + a1
            stdev = sqrt(a2)
            if (3.0_WP*stdev+mean.gt.log(r)) then
               stdev = (log(r) - mean)/3.0_WP
            end if
            vol = 4.0_WP/3.0_WP * pi*r**3
            mean_vol = 4.0_WP/3.0_WP * pi*exp(3.0_WP*mean+4.5_WP*stdev**2)
            if (mean_vol.ge.VFlo) then
               num_part = ceiling(vol/mean_vol)
            else
               num_part = 0
            end if
            if (num_part.lt.2) then
               num_part = 2
            end if
            !print *, "a1 ", a1, " a2 ", a2, " mean ", mean, " stdev ", stdev, " old_vol ", vol, " mean_vol ", mean_vol, " num_part ", num_part
            new_vol = 0.0_WP
            allocate(new_sizes(num_part))
            do n=1,num_part
               call random_number(num1)
               call random_number(num2)
               if (num1 .le. 0.0_WP) num1 = 1e-16_WP
               new_sizes(n) = random_normal(mean, stdev)
               !new_sizes(n) = sqrt(-2.0_WP*log(num1))*cos(2.0_WP*pi*num2)
               !new_sizes(n) = mean + stdev*new_sizes(n)
               new_sizes(n) = exp(new_sizes(n))
               new_vol = new_vol + 4.0_WP/3.0_WP * pi*new_sizes(n)**3
               !print *, "new size ", new_sizes(n)
            end do
            if (new_vol.ge.VFlo) then
               vol_ratio = (vol/new_vol)**(1.0_WP/3.0_WP)
            else
               vol_ratio = 1.0_WP
            end if

            old_pos = this%lp_spray%p(m)%pos
            old_vel = this%lp_spray%p(m)%vel
            open(newunit=iunit,file=trim('spray-all/droplets'),form='formatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim('spray-all/droplets'))
            do n=1,num_part
               new_sizes(n) = new_sizes(n) * vol_ratio
               print *, "This is a secondary drop of diam ", 2.0_WP*new_sizes(n)
               call random_number(num1)
               call random_number(num2)
               call random_number(num3)
               num1 = (num1 - 0.5_WP)*4.0_WP*this%cfg%dx(indx(1))
               num2 = (num2 - 0.5_WP)*4.0_WP*this%cfg%dy(indx(2))
               num3 = (num3 - 0.5_WP)*4.0_WP*this%cfg%dz(indx(3))
               shift = [num1,num2,num3]
               indx2 = this%cfg%get_ijk_global(this%lp_spray%p(m)%pos + shift,indx)
               call random_number(num1)
               call random_number(num2)
               call random_number(num3)
               num1 = (num1 - 0.5_WP)
               num2 = (num2 - 0.5_WP)
               num3 = (num3 - 0.5_WP)
               vel_shift = cross_product(normalize(old_vel),normalize([num1,num2,num3]))*(r/t_bu)
               if (n.eq.1) then
                  this%lp_spray%p(m)%id  =int(4,8)
                  this%lp_spray%p(m)%d   =2.0_WP*new_sizes(n)
                  this%lp_spray%p(m)%pos =old_pos + shift
                  this%lp_spray%p(m)%vel =old_vel + vel_shift
                  this%lp_spray%p(m)%ind =indx2
                  this%lp_spray%p(m)%flag=0
                  this%lp_spray%p(m)%dt  =0.0_WP
                  this%lp_spray%p(m)%Acol=0.0_WP
                  this%lp_spray%p(m)%Tcol=0.0_WP
                  this%lp_spray%p(m)%t   =0.0_WP

                  ! Output diameter, velocity, and position
                  write(iunit,*) this%lp_spray%p(m)%d,this%lp_spray%p(m)%vel(1),this%lp_spray%p(m)%vel(2),this%lp_spray%p(m)%vel(3),&
                  &norm2([this%lp_spray%p(m)%vel(1),this%lp_spray%p(m)%vel(2),this%lp_spray%p(m)%vel(3)]),this%lp_spray%p(m)%pos(1),&
                  &this%lp_spray%p(m)%pos(2),this%lp_spray%p(m)%pos(3),this%lp_spray%p(m)%id,2.0_WP*r
               else
                  ! Increment particle counter
                  this%lp_spray%np_=this%lp_spray%np_+1
                  ! Make room for new drop
                  call this%lp_spray%resize(this%lp_spray%np_)
                  ! Add the drop
                  this%lp_spray%p(this%lp_spray%np_)%id  =int(4,8)
                  this%lp_spray%p(this%lp_spray%np_)%d   =2.0_WP*new_sizes(n)
                  this%lp_spray%p(this%lp_spray%np_)%pos =old_pos + shift
                  this%lp_spray%p(this%lp_spray%np_)%vel =old_vel + vel_shift
                  this%lp_spray%p(this%lp_spray%np_)%ind =indx2
                  this%lp_spray%p(this%lp_spray%np_)%flag=0
                  this%lp_spray%p(this%lp_spray%np_)%dt  =0.0_WP
                  this%lp_spray%p(this%lp_spray%np_)%Acol=0.0_WP
                  this%lp_spray%p(this%lp_spray%np_)%Tcol=0.0_WP
                  this%lp_spray%p(this%lp_spray%np_)%t   =0.0_WP

                  ! Output diameter, velocity, and position
                  write(iunit,*) this%lp_spray%p(this%lp_spray%np_)%d,this%lp_spray%p(this%lp_spray%np_)%vel(1),this%lp_spray%p(this%lp_spray%np_)%vel(2),this%lp_spray%p(this%lp_spray%np_)%vel(3),&
                  &norm2([this%lp_spray%p(this%lp_spray%np_)%vel(1),this%lp_spray%p(this%lp_spray%np_)%vel(2),this%lp_spray%p(this%lp_spray%np_)%vel(3)]),this%lp_spray%p(this%lp_spray%np_)%pos(1),&
                  &this%lp_spray%p(this%lp_spray%np_)%pos(2),this%lp_spray%p(this%lp_spray%np_)%pos(3),this%lp_spray%p(this%lp_spray%np_)%id,2.0_WP*r
               end if
               
               ! Increment monitoring variables
               if (n.gt.1) this%np_drop=this%np_drop+1
               this%lp_spray%np_new=this%lp_spray%np_new+1
               this%lp_spray%vp_new=this%lp_spray%vp_new+4.0_WP/3.0_WP*pi*new_sizes(n)**3
            end do
            close(iunit)

            deallocate(new_sizes)
         end if
      end do

      ! Synchronize particles
      call this%lp_spray%sync()

   end subroutine secondary_break

end module detection_class