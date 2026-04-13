module vfs_data_class
   use precision,      only: WP
   use config_class,   only: config
   use string,         only: str_medium
   use iterator_class, only: iterator
   use irl_fortran_interface
   implicit none
   private
   
   public :: VFhi,VFlo
   real(WP), parameter :: VFlo=1.0e-12_WP  !< Minimum VF value considered
   real(WP), parameter :: VFhi=1.0_WP-VFlo !< Maximum VF value considered
   
   ! Expose type/constructor/methods
   public :: vfs_base,bcond
   
   ! List of known available bcond for this solver
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   
   ! List of available interface reconstructions schemes for VF
   integer, parameter, public :: lvira=1             !< LVIRA scheme
   integer, parameter, public :: elvira=2            !< ELVIRA scheme
   integer, parameter, public :: mof=3               !< MOF scheme
   integer, parameter, public :: wmof=4              !< Wide-MOF scheme
   integer, parameter, public :: r2p=5               !< R2P scheme
   integer, parameter, public :: youngs=6            !< Youngs' scheme
   integer, parameter, public :: lvlset=7            !< Levelset-based scheme
   integer, parameter, public :: plicnet=8           !< PLICnet
   integer, parameter, public :: r2pnet=9            !< R2Pnet
   integer, parameter, public :: jibben=10           !< PPIC-Jibben
   integer, parameter, public :: cylinder=11         !< Cylinder
   integer, parameter, public :: plic_cylinder=12    !< Cylinder
   integer, parameter, public :: r2p_cylinder=13    !< Cylinder
   
   ! List of available interface transport schemes for VF
   integer, parameter, public :: flux=1             !< Flux-based geometric transport
   integer, parameter, public :: flux_storage=2     !< Flux-based geometric transport with storage of detailed face fluxes
   integer, parameter, public :: remap=3            !< Cell-based geometric transport (faster but fluxes are not available)
   integer, parameter, public :: remap_storage=4    !< Cell-based geometric transport with storage of detailed volume moments
   
   ! IRL cutting moment calculation method
   integer, parameter, public :: recursive_simplex=0 !< Recursive simplex cutting
   integer, parameter, public :: half_edge=1         !< Half-edge cutting (default)
   integer, parameter, public :: nonrecurs_simplex=2 !< Non-recursive simplex cutting
   
   ! Default parameters for volume fraction solver
   integer,  parameter, public :: nband=3                                 !< Number of cells around the interfacial cells on which localized work is performed
   integer,  parameter, public :: advect_band=1                           !< How far we do the transport
   integer,  parameter, public :: distance_band=2                         !< How far we build the distance
   integer,  parameter, public :: max_interface_planes=2                  !< Maximum number of interfaces allowed (2 for R2P)
   real(WP), parameter, public :: volume_epsilon_factor =1.0e-15_WP       !< Minimum volume  to consider for computational geometry (normalized by min_meshsize**3)
   real(WP), parameter, public :: surface_epsilon_factor=1.0e-15_WP       !< Minimum surface to consider for computational geometry (normalized by min_meshsize**2)
   real(WP), parameter, public :: iterative_distfind_tol=1.0e-12_WP       !< Tolerance for iterative plane distance finding
   
   !> Bcond shift value
   integer, dimension(3,6), parameter, public :: shift=reshape([+1,0,0,-1,0,0,0,+1,0,0,-1,0,0,0,+1,0,0,-1],shape(shift))
   
   !> Boundary conditions for the volume fraction solver
   type :: bcond
      type(bcond), pointer :: next                        !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'   !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                     !< Bcond type
      integer :: dir                                      !< Bcond direction (1 to 6)
      type(iterator) :: itr                               !< This is the iterator for the bcond
   end type bcond

   ! Define the BASE type with purely data
   type, abstract :: vfs_base
       class(config), pointer :: cfg                       !< This is the config the solver is build for

       ! Volume fraction data
       real(WP), dimension(:,:,:), allocatable :: VF       !< VF array
       real(WP), dimension(:,:,:), allocatable :: VFold    !< VFold array
       
       ! Phase barycenter data
       real(WP), dimension(:,:,:,:), allocatable :: Lbary  !< Liquid barycenter
       real(WP), dimension(:,:,:,:), allocatable :: Gbary  !< Gas barycenter
       
       ! Superficial fluxing velocities
       real(WP), dimension(:,:,:,:), allocatable :: UFl    !< Superficial liquid fluxing velocity
       real(WP), dimension(:,:,:,:), allocatable :: UFg    !< Superficial gas fluxing velocity
       
       ! Subcell phasic volume fields
       real(WP), dimension(:,:,:,:,:,:), allocatable :: Lvol   !< Subcell liquid volume
       real(WP), dimension(:,:,:,:,:,:), allocatable :: Gvol   !< Subcell gas volume
       
       ! Surface density data
       real(WP), dimension(:,:,:), allocatable :: SD       !< Surface density array
       
       ! Distance level set
       real(WP) :: Gclip                                   !< Min/max distance
       real(WP), dimension(:,:,:), allocatable :: G        !< Distance level set array
       
       ! Curvature
       real(WP), dimension(:,:,:), allocatable :: curv     !< Interface mean curvature
       real(WP), dimension(:,:,:,:), allocatable :: curv2p !< Curvature for each interface

       real(WP), dimension(:,:,:), allocatable :: thickness         !< Local thickness of thin region
       real(WP), dimension(:,:,:), allocatable :: thin_sensor       !< Thin structure sensing (=1 is liquid, =2 is gas)
       
       ! Band strategy
       integer, dimension(:,:,:), allocatable :: band      !< Band to localize workload around the interface
       integer, dimension(:,:),   allocatable :: band_map  !< Unstructured band mapping

       ! Masking info for metric modification
       integer, dimension(:,:,:), allocatable :: mask      !< Integer array used for enforcing bconds
       integer, dimension(:,:,:), allocatable :: vmask     !< Integer array used for enforcing bconds - for vertices

             
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_VFS'     !< Solver name (default=UNNAMED_VFS)
      
      ! Boundary condition list
      integer :: nbc                                      !< Number of bcond for our solver
      type(bcond), pointer :: first_bc                    !< List of bcond for our solver
      
      ! Interface handling methods
      integer :: reconstruction_method                    !< Interface reconstruction method
      integer :: transport_method                         !< Interface transport method
      logical :: cons_correct=.true.                      !< Conservative correction (true by default)
      
      ! Flotsam removal parameter
      real(WP) :: flotsam_thld=0.0_WP                     !< Threshold VF parameter for flotsam removal (0.0=off by default)

      ! Parameters for SGS modeling of thin structures
      logical  :: two_planes                              !< Whether we're using a 2-plane reconstruction approach
      logical  :: ppic                                    !< Whether we're using PPIC interfaces
      logical  :: cyl                                     !< Whether we're using cylinder interfaces
      real(WP) :: twoplane_thld1=0.99_WP                  !< Average normal magnitude threshold for r2p to switch from one-plane to two-planes (purely local)
      real(WP) :: twoplane_thld2=0.5_WP                   !< Average normal magnitude threshold above which r2p switches to LVIRA (based on 3x3x3 stencil)
      real(WP) :: thin_thld_dotprod=-0.5_WP               !< Maximum dot product of two interface normals for their respective cells to be considered thin region cells
      real(WP) :: thin_thld_max=0.8_WP                    !< Maximum local thickness to be considered a thin region cell (as a factor of mesh size)
      real(WP) :: thin_thld_min=0.0_WP                    !< Minimum local thickness for thin structure removal (0.0=off, as a factor of mesh size)
      real(WP) :: edge_thld=0.80_WP                       !< Threshold for classifying edges
      
      ! Curvature clipping parameter
      real(WP) :: maxcurv_times_mesh=1.0_WP               !< Clipping parameter for maximum curvature (classically set to 1, but could be larger with r2p since we resolve more)
      
      ! Interface smoothing parameters
      integer  :: smoothing_maxite=0                      !< Maximum number of interface smoothing steps performed after the reconstruction
      real(WP) :: smoothing_maxres=0.0_WP                 !< Maximum residual for interface smoothing - infinity norm, once reached, smoothing stops
      
      ! IRL objects
      type(ByteBuffer_type) :: send_byte_buffer
      type(ByteBuffer_type) :: recv_byte_buffer
      type(ObjServer_SeparatorVariant_type)  :: planar_separator_allocation
      type(ObjServer_PlanarLoc_type)  :: planar_localizer_allocation
      type(ObjServer_LocVariantLink_type) :: localized_separator_link_allocation
      type(ObjServer_LocLink_type)    :: localizer_link_allocation
      type(ObjServer_MixedPolygonBezierSurface_type)    :: interface_mixed_surface_allocation
      type(PlanarLoc_type),        dimension(:,:,:),   allocatable :: localizer
      type(SeparatorVariant_type), dimension(:,:,:),   allocatable :: liquid_gas_interface
      type(LocVariantLink_type),   dimension(:,:,:),   allocatable :: localized_separator_link
      type(ListVM_VMAN_type),      dimension(:,:,:),   allocatable :: triangle_moments_storage
      type(LocLink_type),          dimension(:,:,:),   allocatable :: localizer_link
      type(Poly_type),             dimension(:,:,:,:), allocatable :: interface_polygon
      type(Poly_type),             dimension(:,:,:,:), allocatable :: polyface
      type(SepVM_type),            dimension(:,:,:,:), allocatable :: face_flux    !< Only stored if flux-based transport is used
      type(MixedPolygonBezierSurface_type),  dimension(:,:,:), allocatable :: interface_mixed_surface !< Only used for writing surface (for now!)
      
      ! Monitoring quantities
      real(WP) :: VFmax,VFmin,VFint,SDint                 !< Maximum, minimum, and integral volume fraction and surface density
      real(WP) :: flotsam_error                           !< Integral of flotsam removal error
      real(WP) :: thinstruct_error                        !< Integral of thin structure removal error
      
      integer, dimension(0:nband) :: band_count           !< Number of cells per band value

      ! Additional storage option for geometric transport of auxiliary quantities
      type(TagAccVM_SepVM_type), dimension(:,:,:,:), allocatable :: detailed_face_flux    !< Cell-decomposed face flux geometric data
      type(TagAccVM_SepVM_type), dimension(:,:,:),   allocatable :: detailed_remap        !< Cell-decomposed remapped cell geometric data
      
      ! Old arrays that are needed for the compressible MAST solver
      real(WP), dimension(:,:,:,:), allocatable :: Lbaryold  !< Liquid barycenter
      real(WP), dimension(:,:,:,:), allocatable :: Gbaryold  !< Gas barycenter
      type(SeparatorVariant_type),  dimension(:,:,:), allocatable :: liquid_gas_interfaceold
      type(LocVariantLink_type), dimension(:,:,:), allocatable :: localized_separator_linkold
      type(ObjServer_SeparatorVariant_type)  :: planar_separatorold_allocation
      type(ObjServer_LocVariantLink_type) :: localized_separator_linkold_allocation

   contains
      procedure(sync_interface), deferred :: sync_interface                             !< Synchronize the IRL objects
      procedure(clean_irl_and_band), deferred :: clean_irl_and_band                     !< After a manual change in VF (maybe due to transfer to drops), update IRL and band info
   end type vfs_base

   abstract interface
      subroutine clean_irl_and_band(this)
         import :: vfs_base
         class(vfs_base), intent(inout) :: this
      end subroutine clean_irl_and_band

      subroutine sync_interface(this)
         import :: vfs_base
         class(vfs_base), intent(inout) :: this
      end subroutine sync_interface
   end interface
end module vfs_data_class