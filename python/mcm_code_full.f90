! FFLAGS="-fopenmp -fPIC -Ofast -ffree-line-length-none" f2py-2.7 -c -m mcm_code_full mcm_code_full.f90 wigner3j_sub.f -lgomp
! scinet: FFLAGS="-qopenmp -Ofast -fPIC -xhost" f2py --fcompiler=intelem --noopt -c -m mcm_code_full mcm_code_full.f90 wigner3j_sub.f -liomp5

subroutine calc_mcm(cl,cl_TP,cl_PT,cl_pol,wl_TT,wl_TP,wl_PT,wl_PP, mcm, mcm_Tp, mcm_pT, mcm_pp, mcm_mm)
    implicit none
real(8), intent(in)    :: cl(:),cl_TP(:),cl_PT(:),cl_pol(:),wl_TT(:),wl_TP(:),wl_PT(:),wl_PP(:)
    real(8), intent(inout) :: mcm(:,:),mcm_Tp(:,:),mcm_pT(:,:),mcm_pp(:,:),mcm_mm(:,:)
    real(8), parameter     :: pi = 3.14159265358979323846264d0
    integer :: l1, l2, l3, info, nlmax, lmin, lmax, i
    real(8) :: l1f(2), fac_TT,fac_TP,fac_PT,fac_PP
    real(8) :: thrcof0(2*size(mcm,1)),thrcof1(2*size(mcm,1))     
    nlmax = size(mcm,1)-1

    ! l2,l3 are running over 2,lmax, note that it's tempting to make them run over 0,lmax but the wigner3j associated with polarisation
    ! crash for l2,l3 <0 ( |l2-m2| & |l3-m3| has to be greater or equal to zero)
    ! mcm(l2=2,l3=2)=mcm(1,1)
    ! l1 run in all the non zero wigner 3j it can start at 0 so cl(l=n)=cl(n+1) so cl(l=0)=cl(1)

    !$omp parallel do private(l3,l2,l1,fac_TT,fac_TP,fac_PT,fac_PP,info,l1f,thrcof0,thrcof1,lmin,lmax,i) schedule(dynamic)
    do l1 = 2, nlmax
        write(*,*) l1
        fac_TT=(2*l1+1)/(4*pi)*wl_TT(l1+1)
        fac_TP=(2*l1+1)/(4*pi)*wl_TP(l1+1)
        fac_PT=(2*l1+1)/(4*pi)*wl_PT(l1+1)
        fac_PP=(2*l1+1)/(4*pi)*wl_PP(l1+1)

        do l2 = 2, nlmax

            call drc3jj(dble(l1),dble(l2),0d0,0d0,l1f(1),l1f(2),thrcof0, size(thrcof0),info)
            call drc3jj(dble(l1),dble(l2),-2d0,2d0,l1f(1),l1f(2),thrcof1, size(thrcof1),info)

            lmin=INT(l1f(1))
            lmax=MIN(nlmax,INT(l1f(2)))

            do l3=lmin,lmax
                i   = l3-lmin+1
                mcm(l1-1,l2-1) =mcm(l1-1,l2-1)+ fac_TT*(cl(l3+1)*thrcof0(i)**2d0)
                mcm_Tp(l1-1,l2-1) =mcm_Tp(l1-1,l2-1)+ fac_TP*(cl_TP(l3+1)*thrcof0(i)*thrcof1(i))
                mcm_pT(l1-1,l2-1) =mcm_pT(l1-1,l2-1)+ fac_PT*(cl_PT(l3+1)*thrcof0(i)*thrcof1(i))
                mcm_pp(l1-1,l2-1) =mcm_pp(l1-1,l2-1)+ fac_PP*(cl_pol(l3+1)*thrcof1(i)**2*(1+(-1)**(l1+l2+l3))/2)
                mcm_mm(l1-1,l2-1) =mcm_mm(l1-1,l2-1)+ fac_PP*(cl_pol(l3+1)*thrcof1(i)**2*(1-(-1)**(l1+l2+l3))/2)
    
            end do
        end do
    end do

end subroutine

subroutine bin_mcm(mcm, binLo,binHi, binsize, mbb,scalefac)
    ! Bin the given mode coupling matrix mcm(0:lmax,0:lmax) into
    ! mbb(nbin,nbin) using bins of the given binsize
    implicit none
    real(8), intent(in)    :: mcm(:,:)
    integer, intent(in)    :: binLo(:),binHi(:),binsize(:),scalefac
    real(8), intent(inout) :: mbb(:,:)
    integer :: b1, b2, l1, l2, lmax
    lmax = size(mcm,1)-1
    mbb  = 0
    
    do b2=1,size(mbb,1)
        do b1=1,size(mbb,1)
            do l2=binLo(b2),binHi(b2)
                do l1=binLo(b1),binHi(b1)
                    if (scalefac .eq. 1) then
                        mbb(b1,b2)=mbb(b1,b2) + mcm(l1-1,l2-1)*l2*(l2+1d0)/(l1*(l1+1d0)) !*mcm(l2-1,l3-1)
                    else
                        mbb(b1,b2)=mbb(b1,b2) + mcm(l1-1,l2-1)
                    end if
                end do
            end do
            mbb(b1,b2) = mbb(b1,b2) / binsize(b2)

        end do
    end do
end subroutine

subroutine binning_matrix(mcm, binLo,binHi, binsize, bbl,scalefac)
    implicit none
    real(8), intent(in)    :: mcm(:,:)
    integer(8), intent(in)    :: binLo(:),binHi(:),binsize(:),scalefac
    real(8), intent(inout) :: bbl(:,:)
    integer(8) :: b2, l1, l2,lmax

    lmax = size(mcm,1)-1
    ! mcm is transposed
    ! compute \sum_{l'} M_l'l
    do l1=2,lmax
        do b2=1,size(binLo)
            do l2=binLo(b2),binHi(b2)
                if (scalefac .eq. 1) then
                    bbl(l1-1,b2)=bbl(l1-1,b2)+mcm(l1-1,l2-1)*l2*(l2+1d0)/(l1*(l1+1d0))
                else
                     bbl(l1-1,b2)=bbl(l1-1,b2)+mcm(l1-1,l2-1)
                end if
            end do
            bbl(l1-1,b2)=bbl(l1-1,b2)/(binsize(b2)*1d0)
        end do
    end do

end subroutine



