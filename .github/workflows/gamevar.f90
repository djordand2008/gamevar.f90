!
! gamevar.f90:             calculates the (co)variances of gametic diversity
!
!                        written by Santos, D.J.A  2019
!                        Departament of Animal and Avian Sciences, University of Maryland
!                        daniel_jordan2008@hotmail.com
!
! PARAMETE FILE DEFINITIONS:
! NUMBER OF TRAITS
! interger
! MAP_FILE
! file1
! ALLELE_EFFECTS_FILE
! file2
! PHASED_GENOTYPE_FILE
! file3
! RECOMBINATION_GROUP_FILE
! file4
! GAMETIC_VAR
! T or TRUE
! GAMETIC_(CO)VAR
! T or TRUE
! HOM
! T or TRUE
! CRV
! T or TRUE
! GEBV
! T or TRUE
! OUTPUT_NAME
! name
! OPTION PLINK

program gamevar

implicit none

!-------------------- presentation ---------------------!
character(5):: version= '1.0'
character(13):: date='10 May 2019'

!----------------- parameter file reads ----------------!
integer:: ntraits
character(1)::line_G
character(10),allocatable::group_re(:)
character(22):: parfile,map_imput,allele_imput,geno_imput,re_group_imput,outputfile,reco_imput
logical::C_gamevar=.false.,C_gamecovar=.false.,C_hom=.false.,C_crv=.false.,C_ebv=.false.,C_reco=.false.
logical::C_PLINK=.false.

!----------------- imput variables----------------------!
integer:: imp,io_stat,n_snp_G,n_snp,n_eff,n_anim,n_hap,chr,position,&
i,j,k,t,t1,reco_count,reco_count1,traits_count,traits_count1,&
br,nc_PED,SEX,PED_F,PED_M
integer,allocatable::matrix_count(:),ID_vec(:),group_re_c(:)
real::ID_R,PED_P
real,allocatable::reco(:,:),reco_before(:),traits(:,:)
character(40)::marker,PED_FA  ! changes according to the size
character(200)::outformat=''

!---------------- erro variables -----------------------!
integer::pos_before=-1,chr_before,ID_G
character(len=:), allocatable :: test_p,te_pp
logical::nsort_status=.false.,dchr_status=.false.,nsort_rec_status=.false.,&
colrec_ex_status=.false.,coltrait_ex_status=.false.,n_id_stat=.false.

!------------------ analyses variables ----------------!
integer::re,n, a,cont_sta
integer,allocatable::genotype(:),genobefore(:),genotypeT(:)
real,allocatable::matri_alle_effe(:,:),ebv_results(:),reco_mat_k(:,:),reco_matrix_resu(:,:),&
m_phase(:,:),result_CRV(:),result_HOM(:)
real*10,allocatable::result_GAV(:)

!---------------- output variables --------------------!
integer::nchar,nchar2,ncharebv,nchar2ebv,ncharvar,nchar2var
character(10)::fake ! changes according to the size
character(40)::fake1,fake1_ebv ! changes according to the size
character(80)::fake1_var
character(200)::linhaCRV,linhaEBV,linhaHOM  ! changes according to the size
character(800)::linhaVAR ! changes according to the size

!-------------------------------------------------------!

call print_v
call get_command_argument(1,value=parfile,status=imp)

if(imp/=0) then
    print *,'PARAMETE FILE:'
    read(*,*)parfile
endif

open(unit=2,file=parfile,status='old')

call read_par_file()
call inquire_file()

open(unit=3,file=map_imput,status='old') !chr,marker,position,reco => No head
open(unit=4,file=allele_imput,status='old') ! effects => No head
if (.not. C_PLINK) then
    open(unit=5,file=re_group_imput,status='old')
endif
open(unit=7,file=geno_imput,status='old') ! ID, genotype => No head

n_snp=0
do
    read(3,*,iostat=io_stat) chr !chr,marker,position,reco => No head
    if(io_stat/=0) exit
    n_snp=n_snp+1
enddo
rewind(3)

n_eff=0
do
    read(4,*,iostat=io_stat) ! effects => no head
    if(io_stat/=0) exit
    n_eff=n_eff+1
enddo
rewind(4)

if(n_snp==0 .or. n_eff==0 .or. n_snp/=n_eff) then
	print *
	print*,'*------------------  WARNING -------------------------------*'
	print *
        if(n_snp==0) print*,'    Problem to read the SNP file   '
	if(n_eff==0) print*,'    Problem to read the allele effects file   '
	if(n_snp/=n_eff) print*,'    Numbers of lici with effects and SNPs are different   '
	print *
	print*,'    gamevar.f90 has sttoped    '
	print *
	print*,'*------------------------------------------------------------*'
	print *
endif

allocate(traits(n_snp,ntraits))
traits_count=0
traits_count1=0
do I=1,n_eff
	read(4,*,iostat=io_stat) traits(I,:)
    	if(io_stat/=0) exit
	traits_count1=traits_count1+1
   	if(traits(I,ntraits)==1) then;
        traits_count=traits_count+1
    endif
enddo
if(traits_count==n_eff .or. traits_count1<n_eff) coltrait_ex_status=.TRUE.
rewind(4)

n_anim=0
if (C_PLINK) then
    do
       read(7,*,iostat=io_stat) !group,ID => no header
        if(io_stat/=0) exit
        n_anim=n_anim+1
    enddo
    rewind(7)
else
	do
		read(5,*,iostat=io_stat) !group,ID
		if(io_stat/=0) exit
		n_anim=n_anim+1
	enddo
	rewind(5)
endif

if (C_PLINK) then
    k=1
else
	allocate(group_re(n_anim))
	allocate(group_re_c(n_anim))
	allocate(matrix_count(n_anim))
	allocate(ID_vec(n_anim))
	do I=1,(n_anim)
		read(5,*)group_re(I),ID_vec(I)
		matrix_count(I)=I
	enddo

	k=0
	t=0
	do
		if (sum(matrix_count)==0) exit
		do i=1,n_anim
           if(matrix_count(i)>0) exit
        enddo
        t=t+1
        if(matrix_count(i)>0) then
            do j=i,MAXVAL(matrix_count)
                if (group_re(i)==group_re(j)) then
                    group_re_c(j)=t
                    matrix_count(j)=0
                endif
            enddo
        endif
        k=k+1
	enddo
	deallocate(matrix_count)
	deallocate(group_re)
	
endif

allocate(reco(n_snp,k))
allocate(reco_before(k))
chr_before=chr
reco_before=0
reco_count=0
reco_count1=0

do I=1,n_snp
    if (C_PLINK) then
        read(3,*,iostat=io_stat) chr,marker,reco(I,1),position ! aqui
    else
		read(3,*,iostat=io_stat) chr,marker,position,reco(I,:)
	endif
    if(io_stat/=0) exit
    if(position<pos_before) nsort_status=.TRUE.
    if(chr_before/=chr) dchr_status=.TRUE.
    if(minval(reco(I,:)-reco_before)<0) nsort_rec_status=.TRUE.
    if (position==pos_before) then
        print *
        print*,'*------------------  WARNING ---------------*'
        print *
        print*,'There are markers with same position!'
        print *
        print*,'*-------------------------------------------*'
        print *
    endif
    pos_before=position
    chr_before=chr
    reco_count1=reco_count1+1

    if(reco(I,k)==1) then;
        reco_count=reco_count+1
    endif
    reco_before=reco(I,:)
    enddo
deallocate(reco_before)
if(reco_count==n_snp .or. reco_count1<n_snp) colrec_ex_status=.TRUE.

n_snp_G=0;br=0; nc_PED=0 
do
    read(7,'(a)',advance='no',iostat=io_stat) line_G
    if (is_iostat_eor(io_stat)) then
        exit
    else
        n_snp_G = n_snp_G+1
		if (C_PLINK .and. line_g==' ') then; br=br+1
            if(br==6) nc_PED=n_snp_G-1
        endif
    end if
end do
rewind(7)

if (C_PLINK) then
    n_snp_G=((n_snp_G-nc_PED)/2)/2
else 
	read(7,*,iostat=io_stat) ID_R
	rewind(7)
	n_snp_G=n_snp_G-ceiling(log10(ID_R+1))-1
endif


allocate(matri_alle_effe(n_snp,ntraits))

do I=1,n_eff
read(4,*,iostat=io_stat) matri_alle_effe(I,:)
    if(io_stat/=0) exit
enddo

call print_inf()
call warnings()

n_hap=0
allocate(genotype(n_snp))
allocate(genobefore(n_snp))
if(C_PLINK) allocate(genotypeT(2*n_snp))

allocate(ebv_results(ntraits))
allocate(result_CRV(ntraits))
allocate(result_HOM(ntraits))
if ((C_crv .or. C_gamecovar) .and. .NOT.(C_gamecovar)) then
    allocate(result_GAV(ntraits))
endif

if(C_gamecovar .or. C_crv) then
    allocate(result_GAV(ntraits+(ntraits*ntraits-ntraits)/2))
endif

if (C_gamevar .or. C_gamecovar .or. C_crv) then
    allocate(reco_mat_k(n_snp,n_snp))
    allocate(m_phase(n_snp,n_snp))
    if(k<7) then
        allocate(reco_matrix_resu(n_snp*k,n_snp))
    else
        allocate(reco_matrix_resu(n_snp,n_snp))
    endif
endif

do re=1,k
    call reco_matrix(re,reco_imput)
        reco_matrix_resu(((re-1)*n_snp+1):(re*n_snp),1:n_snp)=reco_mat_k
enddo

if(k<7) then
    deallocate(reco_mat_k)
endif

call head()

if(C_ebv) then
    open (unit=1020,file=trim(outputfile)//'_EBV',status='replace')
    write(1020,*) 'ID '//test_p
endif

if(C_crv) then
    open (unit=1005,file=trim(outputfile)//'_CRV',status='replace')
    write(1005,*) 'ID '//test_p
endif

if(C_hom) then
    open (unit=1010,file=trim(outputfile)//'_HOM',status='replace')
    write(1010,*) 'ID '//test_p
endif

if (C_gamevar .or. C_gamecovar) then
    if((C_gamevar .and. C_gamecovar) .or. C_crv) then
        open (unit=1001,file=trim(outputfile)//'_COVAR',status='replace')
        write(1001,*)'ID '//test_p//te_pp
    endif

    if(C_gamecovar .eqv. .FALSE.) then
        open (unit=1001,file=trim(outputfile)//'_VAR',status='replace')
        write(1001,*)'ID '//test_p
    endif

    if(C_gamevar .eqv. .FALSE.) then
        open (unit=1001,file=trim(outputfile)//'_COVAR',status='replace')
        write(1001,*)'ID '//te_pp
    endif
endif

print*
print*,'Calculating (co)variance of gametic diversity CHR:',chr
print*
i=1
j=1
do
	a=0
    if (C_PLINK) then
        read(7,*,iostat=io_stat) PED_FA,ID_G,PED_F,PED_M,SEX,PED_P,genotypeT(:)
        if(io_stat/=0) exit
        n_hap=n_hap+1
        genotypeT=2*genotypeT-2
        do while (a<(2*n_snp))
            genobefore(a/2+1)=genotypeT(a+2)
            a=a+1
            genotype((a+1)/2)=genotypeT(a)
            a=a+1
        enddo
        write(fake,'(i0)') ID_G
        cont_sta=1
    else
write(outformat,'(a,i0,a,i0,a)') '(i',ceiling(log10(real(ID_vec(i))+1)),',1x,',n_snp,'i1)'
    read(7,outformat,iostat=io_stat) ID_G,genotype(:) ! ID, genotype => No head
        if(io_stat/=0) exit
        n_hap=n_hap+1
        if (Mod((n_hap-1), 2)==0 .and. ID_G/=ID_vec((n_hap/2)+1)) n_id_stat=.TRUE.
        if(Mod((n_hap), 2)==0 .and. ID_G/=ID_vec((n_hap)/2)) n_id_stat=.TRUE.

        if(MOD(j,2) .eq. 0) then
            write(fake,'(i0)') ID_vec(i)
			cont_sta=1
        else
            cont_sta=2
        endif
    endif

    if(cont_sta==1) then
			
        if (C_gamevar .or. C_gamecovar .or. C_crv) then
            do n=1,n_snp
                m_phase(:,n)=(genotype(n)-1)*(genotype-1)
            enddo
			
			if (C_PLINK) then
                m_phase=m_phase*reco_matrix_resu(1:n_snp,1:n_snp)
            else
				if (k<7) then
					m_phase=m_phase*reco_matrix_resu(((group_re_c(I)-1)*n_snp+1):(k*n_snp),1:n_snp)
				else
					call reco_matrix(group_re_c(I),reco_imput)
					m_phase=m_phase*reco_mat_k
				endif
			endif
			do n=1,n_snp
				m_phase(n,n)=0.25
				if(genotype(n)==genobefore(n)) then
					m_phase(n,:)=0
					m_phase(:,n)=0
				endif
            enddo
            if(C_gamecovar) then
                re=ntraits+1
                do t=1,ntraits
                    do t1=1,ntraits
                        if(t==t1) then
                            result_GAV(t)=dot_product(matmul(traits(:,t),m_phase),traits(:,t))
                        endif
                        if(t<t1) then
                            result_GAV(re)=dot_product(matmul(traits(:,t),m_phase),traits(:,t1))
                            re=re+1
                        endif
                    enddo
                enddo
            else
                do t=1,ntraits
                    result_GAV(t)=dot_product(matmul(traits(:,t),m_phase),traits(:,t))
                enddo
            endif

            if (C_gamevar .or. C_gamecovar) then
                ncharvar=0
                ncharvar=len_trim(ADJUSTL(fake))
                linhaVAR(1:(ncharvar+1))=fake
                ncharvar=ncharvar+2

                do t=1,(ntraits+(ntraits*ntraits-ntraits)/2)
                    write(fake1_var,'(F15.5)') result_GAV(t) ! changes according to the size
                    nchar2var=len_trim(ADJUSTL(fake1_var))
                    linhaVAR(ncharvar:(ncharvar+nchar2var))=ADJUSTL(fake1_var)
                    ncharvar=ncharvar+nchar2var+1
                enddo

                write(1001,*) linhaVAR(1:ncharvar-2)
                if(i==1 .or. i==10 .or. i==1000 .or. i==1000000 .or. i==1000000) then
                    print*,'Animal  ..........................................',i
                endif
            endif
        endif

        genobefore=((genobefore+genotype)/2)-1

        if (C_ebv) then
            do t=1,ntraits
                ebv_results(t)=dot_product(genobefore,matri_alle_effe(:,t))
            enddo
            ncharebv=0
            ncharebv=len_trim(ADJUSTL(fake))
            linhaEBV(1:(ncharebv+1))=fake
            ncharebv=ncharebv+2
            do t=1,ntraits
                write(fake1_ebv,'(F15.5)') ebv_results(t) ! changes according to the size
                nchar2ebv=len_trim(ADJUSTL(fake1_ebv))
                linhaEBV(ncharebv:(ncharebv+nchar2ebv))=ADJUSTL(fake1_ebv)
                ncharebv=ncharebv+nchar2ebv+1
            enddo
            write(1020,*) linhaEBV(1:ncharebv-2)
        endif

	    genobefore=genobefore*genobefore

        if (C_hom) then
            do t=1,ntraits
                result_HOM(t)=dot_product(genobefore,(matri_alle_effe(:,t)*matri_alle_effe(:,t)))
            enddo
            nchar=0
            nchar=len_trim(ADJUSTL(fake))
            linhaHOM(1:(nchar+1))=fake
            nchar=nchar+2
            do t=1,ntraits
                write(fake1,'(F10.2)') result_HOM(t) ! changes according to the size
                nchar2=len_trim(ADJUSTL(fake1))
                linhaHOM(nchar:(nchar+nchar2))=ADJUSTL(fake1)
                nchar=nchar+nchar2+1
            enddo
            write(1010,*)linhaHOM(1:nchar-2)
        endif

		if (C_crv) then
            if(C_hom) then
				do t=1,ntraits
					result_CRV(t)=sqrt(result_GAV(t))/(sqrt(0.5*(result_HOM(t))+result_GAV(t)))
				enddo
			else
				do t=1,ntraits
                    result_CRV(t)=dot_product(genobefore,(matri_alle_effe(:,t)*matri_alle_effe(:,t)))
                    result_CRV(t)=sqrt(result_GAV(t))/(sqrt(0.5*(result_CRV(t))+result_GAV(t)))
                enddo
			endif
            nchar=0
            nchar=len_trim(ADJUSTL(fake))
            linhaCRV(1:(nchar+1))=fake
            nchar=nchar+2
            do t=1,ntraits
                write(fake1,'(F10.2)') result_CRV(t) ! changes according to the size
                nchar2=len_trim(ADJUSTL(fake1))
                linhaCRV(nchar:(nchar+nchar2))=ADJUSTL(fake1)
                nchar=nchar+nchar2+1
            enddo
            write(1005,*)linhaCRV(1:nchar-2)
        endif
            i=i+1
    endif
    if(.not. C_PLINK) genobefore=genotype
    if(.not. C_PLINK) j=j+1
enddo
rewind(7)
if(C_PLINK) then
    deallocate(genotypeT)
    n_hap=2*n_hap
endif
deallocate(reco_matrix_resu)
print *
print *,'Number of Haplotypes analysed:....................',n_hap

print *
print *, 'Outputs:'
print *, '------------------------------------------------------'

if (C_gamevar .or. C_gamecovar) then
    if (C_gamecovar .eqv. .false.) then
        print'(a)','   '//trim(outputfile)//'_VAR has been created!'
    else
        print'(a)','   '//trim(outputfile)//'_COVAR has been created!'
    endif
endif

if(C_hom) then
    print'(a)','   '//trim(outputfile)//'_HOM has been created!'
endif

if(C_crv) then
    print'(a)','   '//trim(outputfile)//'_CRV has been created!'
endif

if(C_ebv) then
    print'(a)','   '//trim(outputfile)//'_EBV has been created!'
endif

print *, '------------------------------------------------------'

print *
print *
call printime()
print *
print *
print *, 'The solutions have been calculated!!! Have a nice day !!!      ;)'
print *

CONTAINS

!-------------------------------print_v-------------------------------------!
subroutine print_v()
print *
print *,'  *---------------------------------------------------------------------------*'
print *,'  *----------------------   gamevar.f90 ','Version: ', version, '------------------------*'
print *,'  *---------------------------------------------------------------------------*'
print *,'  * Authorship: Daniel Jordan de Abreu Santos                                 *'
print *,'  * Algorithm basead on Santos, D.J.A et al., Journal of Dairy Sciences, 2019 *'
print *,'  * Creation data and location: ',date,'- UMD, USA','                       *'
print *,'  *---------------------------------------------------------------------------*'
print *
call printime()
print *
print *
end subroutine

!-------------------------------printime-------------------------------------!
subroutine printime

implicit none

integer  DATE_TIME (8)
character(LEN = 12) REAL_CLOCK (3)

call DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), &
REAL_CLOCK (3), DATE_TIME)

print *,'Current Data and Time: ',real_clock(1)(7:8),'/',real_clock(1)(5:6),&
'/',real_clock(1)(1:4),' ', &
real_clock(2)(1:2),':',real_clock(2)(3:4),':',real_clock(2)(5:6)
end subroutine

!-------------------------------print_inf-------------------------------------!
subroutine print_inf()

implicit none

print *
print *
print *,'DESCRIPTION:'
print *
print *,'Parameter File Readed:      ',parfile
print *,'Number of SNPs:          ',n_snp
print *,'Number of Allele Effects:',n_eff
print *,'Number of Genotypes:     ',n_snp_G
print *,'Number of Reco Groups:   ',k
print *,'Number of Traits:        ',ntraits
print *,'Number of Animals:       ',n_anim
print *
end subroutine

!------------------------------read_par_file--------------------------------------!
subroutine read_par_file()

implicit none

integer :: out_test,io_statP 
character(4):: stat_var,stat_covar,stat_hom,stat_crv,stat_ebv
character(6)::OPTI
character(12)::stat_PLINK

print *
print *, 'OPTIONS:'
print *
read(2,*)
read(2,*) ntraits;print '(A,I4)',' Number of Traits: ',ntraits
read(2,*)
read(2,*) map_imput;print *,'Map File: ',map_imput
read(2,*)
read(2,*) allele_imput;print *,'Allele Effects File: ',allele_imput
read(2,*)
read(2,*) geno_imput;print *,'Phased Genotype File: ',geno_imput
read(2,*)
read(2,*)re_group_imput
if (re_group_imput=='GENETIC_DISTANCE_UNIT') then; re_group_imput='NO'; BACKSPACE(2);endif
print *,'Recombination Group File: ',re_group_imput
read(2,*)
read(2,*) reco_imput; if(reco_imput=='morgans') then
        C_reco=.true.
        print *,'Genetic Distance Unit: Morgans'
    else if (reco_imput=='reco') then
        C_reco=.true.
        print *,'Genetic Distance Unit: Recombination rate'
    else
        C_reco=.false.
        print*, 'Genetic Distance Unit: No Known Unit'
    endif
read(2,*)
print *
read(2,*)stat_var; if (stat_var=='T' .OR. stat_var=='TRUE') then;C_gamevar=.true.;print *, 'Estimate Gametic Variance = TRUE'; endif
read(2,*)
read(2,*)stat_covar; if (stat_covar=='T' .OR. stat_covar=='TRUE') then
C_gamecovar=.true.;print *, 'Estimate Gametic (Co)Variance = TRUE'; endif
read(2,*)
read(2,*)stat_hom; if (stat_hom=='T' .OR. stat_hom=='TRUE') then;C_hom=.true.;print *, 'Estimate HOM = TRUE'; endif
read(2,*)
read(2,*)stat_crv; if (stat_crv=='T' .OR. stat_crv=='TRUE') then;C_crv=.true.;print *, 'Estimate CRV = TRUE'; endif
read(2,*)
read(2,*)stat_ebv; if (stat_ebv=='T' .OR. stat_ebv=='TRUE') then;C_ebv=.true.;print *, 'Estimate EBV = TRUE'; endif
read(2,*)
print *
read(2,*,IOSTAT=out_test) outputfile; if (out_test==-1) then
    outputfile='youforgotitFDP';endif
print*,'output file: ',outputfile
read(2,*,iostat=io_statP)OPTI,stat_PLINK
print *
if(io_statP==0 .and. OPTI== "OPTION" .and. stat_PLINK=='PLINK') then
    print *
    print *, OPTI,' ',stat_PLINK
    print *
    C_PLINK=.true.; print *, 'Reading PLINK Format ...'; 
endif
if(re_group_imput=='NO' .and. io_statP/=0) then
    print *
    print *,' The recombination group file must to be provided'
print *
endif
end subroutine

!-------------------------------inquire files-------------------------------------!
subroutine inquire_file()

implicit none

logical:: exi_map,exi_alle,geno_exi,group_exi
inquire (file=map_imput,exist=exi_map)
inquire (file=allele_imput,exist=exi_alle)
inquire (file=geno_imput,exist=geno_exi)
inquire (file=re_group_imput,exist=group_exi)
if (C_PLINK) group_exi=.TRUE.

if(.NOT. (exi_map) .OR. .NOT. (exi_alle) .OR. .NOT. (geno_exi) .OR. .NOT. (group_exi)) then
    print *
    print*,'*------------------  WARNING ---------------*'
    print *

    if (.NOT.(exi_map)) then
        print*, map_imput,'is missing'
    endif

    if (.NOT.(exi_alle)) then
        print*, allele_imput,'is missing'
    endif

    if (.NOT.(geno_exi)) then
        print*, geno_imput,'is missing'
    endif

    if (.NOT.(group_exi)) then
        print*, re_group_imput,'is missing'
    endif

    print *
    print*,'    gamevar.f90 has sttoped   '
    print *
    print*,'*--------------------------------------------*'
    print *

    stop

endif

end subroutine

!-------------------------------warnings-------------------------------------!
subroutine warnings()

implicit none

logical::wa_nsnp=.FALSE.,wa_nani=.FALSE.,wa_geno=.FALSE.

if(n_snp/=n_eff) wa_nsnp=.TRUE.
if(n_anim/=(n_hap/2)) wa_nani=.TRUE.
if(n_snp_G/=n_snp) wa_geno=.TRUE.

if(wa_nsnp .OR. wa_nsnp .OR. nsort_status .OR. dchr_status .OR. nsort_rec_status &
.OR. colrec_ex_status .OR. coltrait_ex_status .OR. n_id_stat .OR. wa_geno) then
    print *
    print*,'*------------------  WARNING ---------------------------------*'
    print *
    if(wa_nsnp) then
        print *,'Number of SNP and effects are diferent'
        print *,map_imput,'=',n_snp,' and ',allele_imput,'=',n_eff
        print *
    end if
    if(wa_nani) then
        print *
        print *,'Number of individuals with haplotypes are diferent'
        print *,re_group_imput,'=',n_anim,' and ',geno_imput,'=',n_hap/2
        print *
    end if

    if (nsort_status) then
        print *
        print*, 'Markers are not properly sorted by bp position!'
        print *
    endif

    if (dchr_status) then
        print *
        print*, 'There are more than one chromosome in:',map_imput
        print *
    endif

    if (nsort_rec_status) then
        print *
        print*,'Markers are not properly sorted by recombination rate!'
        print *
    endif

    if (colrec_ex_status) then
        print *
        print*,'Number of recombination columns exceeds the number of recombination groups!'
        print *
    endif

    if (coltrait_ex_status) then
        print *
        print*,'Parameter number of traits exceeds the number of trait columns!'
        print *
    endif

    if (n_id_stat) then
        print *
        print*,'Expected IDs orders from     ',geno_imput,'     and      ',re_group_imput,' do not match!'
        print *
    endif

    if (wa_geno) then
        print *
        print*,'Number of SNPs and genotypes are diferent!'
        print *,map_imput,'=',n_snp,' and ',geno_imput,'=',n_snp_G
        print *
    endif

    print*,'    gamevar.f90 has sttoped   '
    print *
    print*,'*------------------------------------------------------------*'
    print *
    stop
end if

if ((C_crv .or. C_gamevar .or. C_gamecovar) .and. .not. C_reco) then
    print *
    print*,'*------------------  WARNING --------------------------------*'
    print *
    print *
    print*,'   Genetic Unit is required for Variance of Gametic Diversity or CRV Calculation!'
    print *
    print*,'    gamevar.f90 has sttoped   '
    print *
    print*,'*------------------------------------------------------------*'
    print *
    stop
endif

end subroutine

!-------------------------------head-------------------------------------!
subroutine head()

implicit none

integer::fi
real::ir,jr
character(8)::nada,nadaj
character::traits_covec((ntraits*ntraits-ntraits)/2)

fi=0
do i=1,ntraits
    ir=i
    fi=fi+ceiling(log10(ir+1))
enddo

allocate(character(len=7*ntraits+fi)::test_p)
fi=0
do i=1,ntraits
    write(nada,'(i0)') i
    ir=i
    test_p((fi+1):(fi+7+ceiling(log10(ir+1))))='TRAIT_'//nada
    fi=fi+7+ceiling(log10(ir+1))
enddo

if (C_gamecovar) then
    fi=0
    do i=1,(ntraits-1)
        do j=(i+1),ntraits
            if (i<j) then
                ir=i
                jr=j
                fi=fi+ceiling(log10(ir+1))+ceiling(log10(jr+1))
            endif
        enddo
    enddo

    allocate(character(len=((ntraits*ntraits-ntraits)/2)*8+fi)::te_pp)
    fi=0
    do i=1,(ntraits-1)
        do j=(i+1),ntraits
            if (i<j) then
                write(nada,'(i0)') i
                write(nadaj,'(i0)') j
                ir=i
                jr=j
te_pp((fi+1):(fi+7+ceiling(log10(ir+1))))='TRAIT_'//nada
fi=fi+7+ceiling(log10(ir+1))
te_pp((fi):(fi+ceiling(log10(jr+1))))='_'//nadaj
fi=fi+ceiling(log10(jr+1))
  			if((i+j)<(2*ntraits-1)) then
te_pp((fi+1):(fi+2))=' '
fi=fi+1
		endif
            endif
        enddo
    enddo
endif

end subroutine head

!-------------------------------reco_matrix-------------------------------------!
subroutine reco_matrix(qlq,reco_x)

implicit none

integer::d1,d2,qlq,deno=0
character(6)::reco_x

if (reco_x=='reco') deno=2
if (reco_x=='morgans') deno=200
!deno=2

if (k<7) print*,"Calculating recombination matrix ........... group",qlq

do d1=1,n_snp
    do d2=1,n_snp
        if(d1/=d2) then
            if ((sqrt((reco(d1,qlq)-reco(d2,qlq))**2))>=0.50) then
                reco_mat_k(d1,d2)=0
            else
                reco_mat_k(d1,d2)=(-(sqrt((reco(d1,qlq)-reco(d2,qlq))**2)/deno))+0.25
            endif
        endif
    enddo
enddo

end subroutine reco_matrix

end program gamevar


