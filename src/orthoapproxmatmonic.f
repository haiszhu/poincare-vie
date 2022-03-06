      subroutine orthoapproxmatmonic(nsuccess)
          
      implicit real *8 (a-h,o-z)
      real *8 uvs(2,10000),wts(10000)
      real *8, allocatable :: x(:),y(:),w(:),basis_jk(:),basiskl(:)
      real *8, allocatable :: a_diag(:),vals2coefs(:,:)
      real *8, allocatable :: coefs2coefs(:,:)
      real *8, allocatable :: o_basis(:,:),p_coefs(:,:)
      real *8 nrm_kl ! numerator
      integer aorder,idx_jk,idx_jk_new
      integer add_power(2),power_x_new,power_y_new
      integer, allocatable :: o_basis_tile(:,:), tile_to_idx(:)
      integer, allocatable :: power_x(:), power_y(:)
      integer, allocatable :: xpower(:), ypower(:)

      nsuccess = 1

      aorder = 8
      norder = 8
      npols = (norder+1)*(norder+2)/2

      allocate(x(npols),y(npols),w(npols),basis_jk(npols))
      allocate(o_basis(aorder*(aorder+1)/2,npols))
      allocate(p_coefs(aorder*(aorder+1)/2,aorder*(aorder+1)/2))
      allocate(o_basis_tile(aorder,aorder))
      allocate(tile_to_idx(aorder*(aorder+1)/2))
      allocate(a_diag(aorder*(aorder+1)/2))
      allocate(vals2coefs(aorder*(aorder+1)/2,npols))
      allocate(coefs2coefs(aorder*(aorder+1)/2,aorder*(aorder+1)/2))
      allocate(power_x(aorder**2),power_y(aorder**2))
      allocate(xpower(aorder*(aorder+1)/2),ypower(aorder*(aorder+1)/2))


      call get_vioreanu_nodes(norder,npols,uvs)
      call get_vioreanu_wts(norder,npols,wts)
      x = uvs(1,1:npols)
      y = uvs(2,1:npols)
      w = wts(1:npols)
      o_basis = 1.0
      p_coefs = 0.0
      do i = 1,aorder*(aorder+1)/2
          p_coefs(i,i) = 1.0
      end do
      o_basis_tile(1,1) = 1
      tile_to_idx(1) = 1
      do j = 2,aorder
          do k = 1,j
              if (k .lt. j) then
                  basis_jk = x*o_basis(j*(j-1)/2+k-j+1,:)
              else
                  basis_jk = y*o_basis(j*(j-1)/2,:)
              end if
              o_basis(j*(j-1)/2+k,:) = basis_jk
              o_basis_tile(k,j-k+1) = j*(j-1)/2+k;
              tile_to_idx(j*(j-1)/2+k) = (j-k)*aorder+k;
          end do
          do k = j*(j-1)/2+1,j*(j+1)/2 
              do l = 1,k-1
                  nrm_kl = sum(o_basis(k,:)*o_basis(l,:)*w)
                  p_coefs(l,k) = -nrm_kl/sum(o_basis(l,:)**2.*w)
              end do
              basiskl = matmul(transpose(o_basis(1:k,:)),p_coefs(1:k,k))
              o_basis(k,:) = basiskl
          end do
      end do

      a_diag = matmul(o_basis**2,w)
      do i = 1,aorder*(aorder+1)/2
        vals2coefs(i,:) = o_basis(i,:)*w/a_diag(i)
      end do

      coefs2coefs = 0.0
      do i = 1,aorder*(aorder+1)/2
          coefs2coefs(i,i) = 1.0
      end do
      power_x=reshape(spread((/(i,i=0,aorder-1)/),1,aorder),[aorder**2])
      power_y=reshape(spread((/(i,i=0,aorder-1)/),2,aorder),[aorder**2])
      do j = 1,aorder
          do k = j*(j-1)/2+1,j*(j+1)/2 
              if (k .lt. j*(j+1)/2) then
                  idx_jk = k-j+1
                  add_power = (/1,0/)
              else
                  idx_jk = k-j
                  add_power = (/0,1/)
              end if
              do l = 1,idx_jk
                  power_x_new = power_x(tile_to_idx(l))+add_power(1)
                  power_y_new = power_y(tile_to_idx(l))+add_power(2)
                  idx_jk_new = o_basis_tile(power_y_new+1,power_x_new+1)
                  coefs2coefs(idx_jk_new,k) = coefs2coefs(l,idx_jk)
              end do
          coefs2coefs(1:k,k)=matmul(coefs2coefs(1:k,1:k),p_coefs(1:k,k))
          end do
      end do

      ! eval
      xpower = power_x(tile_to_idx); ypower = power_y(tile_to_idx);

      ! 
      open(2, file = '../orthomat.txt') 
      do i = 1,aorder*(aorder+1)/2
          write(2,*) coefs2coefs(i,:)
      end do  
      close(2) 

      end

      
