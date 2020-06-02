! Fortran library for investigating structure and properties of liquid water
! Works for both bulk water and nearby solutes
! SpherePoints and SphereSurfaceAreas are copied from sim.geom.py
! with minor modifications made.


! Same as in geom.py
! Finds centroid of given atom set
subroutine Centroid(Pos, Ret, N, Dim)
    implicit none
    integer, intent(in) :: N, Dim
    real(8), dimension(N,Dim), intent(in) :: Pos
    real(8), dimension(Dim), intent(out) :: Ret
    Ret = sum(Pos, 1) / real(N)
end subroutine

! Taken from geom.py
subroutine crossProd3(r1, r2, Ret, Dim)
    integer, intent(in) :: Dim
    real(8), dimension(Dim), intent(in) :: r1, r2
    real(8), dimension(Dim), intent(out) :: Ret
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    Ret(1) = r1(2)*r2(3) - r1(3)*r2(2)
    Ret(2) = r1(3)*r2(1) - r1(1)*r2(3)
    Ret(3) = r1(1)*r2(2) - r1(2)*r2(1)
end subroutine

! Same as in geom.py
subroutine reimage(Pos, RefPos, BoxL, ReimagedPos, NPos, Dim)
    implicit none
    integer, intent(in) :: NPos, Dim
    real(8), dimension(NPos, Dim), intent(in) :: Pos
    real(8), dimension(Dim), intent(in) :: RefPos
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(NPos, Dim), intent(out) :: ReimagedPos
    integer :: i
    real(8), dimension(Dim) :: distvec, iBoxL
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    do i = 1, NPos
        distvec = Pos(i,:) - RefPos
        distvec = distvec - BoxL * anint(distvec * iBoxL)
        ReimagedPos(i,:) = RefPos + distvec
    enddo            
end subroutine

! Same as in proteinlib.f90
real(8) function RgWeights(Pos, Weights, NAtom, Dim)
  implicit none
  integer, intent(in) :: NAtom, Dim
  real(8), dimension(NAtom,Dim), intent(in) :: Pos
  real(8), dimension(NAtom), intent(in) :: Weights
  real(8), dimension(Dim) :: Center
  integer :: i
  Center = sum(Pos, 1) / real(NAtom)
  RgWeights = 0.
  do i = 1, NAtom
    RgWeights = RgWeights + Weights(i) * sum((Pos(i,:) - Center)**2)
  enddo
  RgWeights = RgWeights / sum(Weights)
  RgWeights = sqrt(RgWeights)
end function

! Same as in geom.py
! Places points on a sphere to later define SASA
subroutine SpherePoints(N, Points)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N,3), intent(out) :: Points
    real(8) :: off, y, phi, r
    real(8) :: inc 
    real(8), parameter :: pi = 3.1415926535897931D0 
    integer :: k
    inc = pi * (3. - sqrt(5.))
    Points = 0.
    off = 2. / real(N)
    do k = 1, N
        y = real(k-1) * off - 1. + (off * 0.5)
        r = sqrt(max(1. - y*y, 0.))
        phi = real(k-1) * inc
        Points(k,1) = cos(phi)*r
        Points(k,2) = y
        Points(k,3) = sin(phi)*r
    enddo
end subroutine

! Modified to also return Exposed, so can find which points (then atoms) are exposed
subroutine SphereSurfaceAreas(Pos, Radii, Points, nExp, BoxL, Areas, Exposed, NSphere, NPoints, Dim)
    implicit none
    integer, intent(in) :: Dim
    real(8), dimension(NSphere, Dim), intent(in) :: Pos
    real(8), dimension(NSphere), intent(in) :: Radii
    real(8), dimension(NPoints, Dim), intent(in) :: Points
    integer, intent(in) :: nExp
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(NSphere), intent(out) :: Areas
    logical, dimension(NSphere), intent(out) :: Exposed
    real(8), parameter :: pi = 3.141592653589D0 
    integer, intent(in) :: NSphere, NPoints
    integer :: i, j, k
    real(8), dimension(NPoints,Dim) :: ThisPoints
    real(8) :: AreaPerPoint
    real(8), dimension(NSphere) :: RadiiSq
    real(8), dimension(Dim) :: iPos, jPos
    real(8), dimension(Dim) :: distvec, iBoxL
    logical, dimension(NPoints) :: tempExposed
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0) 
    Areas = 0.
    Exposed = .false.
    RadiiSq = Radii*Radii
    do i = 1, NSphere
        iPos = Pos(i,:)
        AreaPerPoint = 4.*pi*Radii(i)**2 / real(NPoints)
        tempExposed = .true.
        do k = 1, NPoints
            ThisPoints(k,:) = Points(k,:) * Radii(i) + iPos
        enddo
        do j = 1, NSphere
            if (i == j) cycle
            jPos = Pos(j,:)
            distvec = jPos - iPos 
            distvec = distvec - BoxL * anint(distvec * iBoxL)
            jPos = iPos + distvec
            if (.not. any(tempExposed)) exit
            !first check if spheres are far from each other
            if (sum((jPos-iPos)**2) > (Radii(i) + Radii(j))**2) cycle
            do k = 1, NPoints
                if (.not. tempExposed(k)) cycle
                if (sum((ThisPoints(k,:) - jPos)**2) < RadiiSq(j)) tempExposed(k) = .false.
            enddo
        enddo
        Areas(i) = AreaPerPoint * real(count(tempExposed))
        if (count(tempExposed) >= nExp) Exposed(i) = .true.
    enddo
end subroutine

! Same as in geom.py
subroutine SphereVolumes(Pos, Radii, dx, Volumes, NSphere, Dim)
    implicit none
    integer, intent(in) :: Dim
    real(8), dimension(NSphere, Dim), intent(in) :: Pos
    real(8), dimension(NSphere), intent(in) :: Radii
    real(8),  intent(in) :: dx
    real(8), dimension(NSphere), intent(out) :: Volumes
    integer, intent(in) :: NSphere
    real(8), dimension(NSphere) :: RadiiSq
    real(8) :: minDistSq, DistSq, dV
    integer :: i,j
    real(8), dimension(Dim) :: Pos2, minPos, maxPos
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    RadiiSq = Radii*Radii
    Volumes = 0.
    dV = dx*dx*dx
    minPos = (/minval(Pos(:,1) - Radii), minval(Pos(:,2) - Radii), minval(Pos(:,3) - Radii)/)
    maxPos = (/maxval(Pos(:,1) + Radii), maxval(Pos(:,2) + Radii), maxval(Pos(:,3) + Radii)/)
    maxPos = maxPos + dx * 0.5    
    !first do a coarse grid check to see which spheres are where
    Pos2 = minPos
    do while (all(Pos2 < maxPos))
        j = 0
        minDistSq = huge(1.d0)
        do i = 1, NSphere
            DistSq = sum((Pos(i,:) - Pos2)**2)
            if (DistSq < minDistSq .and. DistSq < RadiiSq(i)) then
                minDistSq = DistSq
                j = i
            endif
        enddo
        if (j > 0) Volumes(j) = Volumes(j) + dV
        Pos2(1) = Pos2(1) + dx
        do i = 1, 2
            if (Pos2(i) >= maxPos(i)) then
                Pos2(i) = minPos(i)
                Pos2(i+1) = Pos2(i+1) + dx
            endif
        enddo
    enddo   
end subroutine

! Calculates average g(r) for atoms in Pos2 to atoms in Pos1
! Will work best for protein if use SURFACE atoms as Pos1, not all
! Assumes a bulk density of atoms in Pos2 provided... i.e. system can be inhomogeneous
! rather than assuming homogeneously distributed atoms in Pos2 over whole volume of box
subroutine RadialDist(Pos1, Pos2, binwidth, totbins, BulkDens, BoxL, rdf, NPos1, NPos2)
    implicit none
    real(8), dimension(NPos1, 3), intent(in) :: Pos1
    real(8), dimension(NPos2, 3), intent(in) :: Pos2
    real(8), intent(in) :: binwidth
    integer, intent(in) :: totbins
    ! Total number of bins, so max bin value is totbins*binwidth; if larger ignored
    real(8), intent(in) :: BulkDens
    real(8), dimension(3), intent(in) :: BoxL
    real(8), dimension(totbins), intent(out) :: rdf
    integer, intent(in) :: NPos1, NPos2
    real(8), parameter :: pi = 3.141592653589D0
    real(8), dimension(totbins) :: counts
    real(8), dimension(3) :: iPos, jPos, iBoxL, distVec
    real(8) :: dist
    integer :: nbin, i, j, k
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    counts = 0.
    do i = 1, NPos2
        iPos = Pos2(i,:)
        do j = 1, NPos1
            jPos = Pos1(j,:)
            distVec = jPos - iPos
            distVec = distVec - BoxL * anint(distVec * iBoxL)
            dist = sqrt(sum(distVec**2))
            ! Place in correct bin
            nbin = ceiling(dist / binwidth) ! Bin i has all distances less than i but greater than i - 1
            if (nbin <= totbins) then
                counts(nbin) = counts(nbin) + 1.
            endif
            ! If beyond max of totbins*binwidth, just don't include
        enddo
    enddo
    ! Normalize counts correctly
    do k = 1, totbins
        rdf(k) = counts(k) / (NPos1 * BulkDens * (4./3.) * pi * (binwidth**3) * ((k**3) - ((k-1)**3)))
    enddo
    ! Gives g(r) for single snapshot... for whole trajectory, do for each frame and average
end subroutine

subroutine RadialDistSame(Pos, binwidth, totbins, BulkDens, BoxL, rdf, NPos)
    implicit none
    real(8), dimension(NPos, 3), intent(in) :: Pos
    real(8), intent(in) :: binwidth
    integer, intent(in) :: totbins
    ! Total number of bins, so max bin value is totbins*binwidth; if larger ignored
    real(8), intent(in) :: BulkDens
    real(8), dimension(3), intent(in) :: BoxL
    real(8), dimension(totbins), intent(out) :: rdf
    integer, intent(in) :: NPos
    real(8), parameter :: pi = 3.141592653589D0
    real(8), dimension(totbins) :: counts
    real(8), dimension(3) :: iPos, jPos, iBoxL, distVec
    real(8) :: dist
    integer :: nbin, i, j, k
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    counts = 0.
    do i = 1, NPos
        iPos = Pos(i,:)
        do j = i+1, NPos
            jPos = Pos(j,:)
            distVec = jPos - iPos
            distVec = distVec - BoxL * anint(distVec * iBoxL)
            dist = sqrt(sum(distVec**2))
            ! Place in correct bin
            nbin = ceiling(dist / binwidth) ! Bin i has all distances less than i but greater than i - 1
            if (nbin <= totbins) then
                counts(nbin) = counts(nbin) + 1.
            endif
            ! If beyond max of totbins*binwidth, just don't include
        enddo
    enddo
    ! Normalize counts correctly
    do k = 1, totbins
        rdf(k) = counts(k) / (NPos * BulkDens * (4./3.) * pi * (binwidth**3) * ((k**3) - ((k-1)**3)))
    enddo
    ! Gives g(r) for single snapshot... for whole trajectory, do for each frame and average
end subroutine

! Calculate the end-to-end distance of a set of position vectors for a particle 1 and 2.
subroutine Ree(Pos1,Pos2,BoxL,dist,Npos1,Npos2,Dim)
    implicit none
    integer, intent(in) :: NPos1, NPos2, Dim
    real(8), dimension(Npos1, Dim), intent(in) :: Pos1
    real(8), dimension(Npos2, Dim), intent(in) :: Pos2
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(Npos1), intent(out) :: dist
    real(8), dimension(Dim) :: iPos1, iPos2, iBoxL, distVec, distVec_re
    real(8), dimension(Npos1) :: dist_re
    integer :: i
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    do i = 1, NPos1
       iPos1 = Pos1(i,:)
       iPos2 = Pos2(i,:)
       distVec = iPos2-iPos1
       distVec = distVec - BoxL * anint(distVec * iBoxL)
       ! test code do another correction for min image convention
       !distVec = distVec - BoxL * anint(distVec * iBoxL)
       dist(i) = sqrt(sum(distVec**2))
    enddo
end subroutine

! Calculate the end-to-end vector of a set of position vectors for a particle 1 and 2.
subroutine Reevec(Pos1,Pos2,BoxL,dist,Npos1,Npos2,Dim)
    implicit none
    integer, intent(in) :: NPos1, NPos2, Dim
    real(8), dimension(Npos1, Dim), intent(in) :: Pos1
    real(8), dimension(Npos2, Dim), intent(in) :: Pos2
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(Npos1, Dim), intent(out) :: dist
    real(8), dimension(Dim) :: iPos1, iPos2, distVec, iBoxL
    integer :: i
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    do i = 1, NPos1
       iPos1 = Pos1(i,:)
       iPos2 = Pos2(i,:)
       distVec = iPos2-iPos1
       distVec = distVec - BoxL * anint(distVec * iBoxL)
       dist(i,:) = distVec
    enddo
end subroutine

! Write more general function for finding probability distribution of distances in any dimension, with reimaging
! If Pos1 and Pos2 are the same, then divide the number of counts by 2
! If the distance is zero, it's not counted
subroutine PairDistanceHistogram(Pos1, Pos2, binwidth, totbins, BoxL, hist, NPos1, NPos2, Dim)
    implicit none
    integer, intent(in) :: NPos1, NPos2, Dim
    real(8), dimension(NPos1, Dim), intent(in) :: Pos1
    real(8), dimension(NPos2, Dim), intent(in) :: Pos2
    real(8), intent(in) :: binwidth
    integer, intent(in) :: totbins
    ! Total number of bins, so max bin value is totbins*binwidth; if larger ignored
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(totbins), intent(out) :: hist
    real(8), dimension(Dim) :: iPos, jPos, iBoxL, distVec
    real(8) :: dist
    integer :: nbin, i, j
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    hist = 0.
    do i = 1, NPos1
        iPos = Pos1(i,:)
        do j = 1, NPos2
            jPos = Pos2(j,:)
            distVec = jPos - iPos
            distVec = distVec - BoxL * anint(distVec * iBoxL)
            dist = sqrt(sum(distVec**2))
            if (dist == 0.d0) cycle ! Skip it if it's the same position
            ! Place in correct bin
            nbin = ceiling(dist / binwidth) ! Bin i has all distances less than i but greater than i - 1
            if (nbin <= totbins) then
                hist(nbin) = hist(nbin) + 1.
            endif
            ! If beyond max of totbins*binwidth, just don't include
        enddo
    enddo
end subroutine            
