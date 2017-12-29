  function M = J(x)
      n=numel(x);
      M=sparse(n, n);
      M(1,1)=2*x(1);
      
%      M(1,2)=-1;
%      M(2,1)=1;
%      M(2,2)=-1;
      
      
      for i=2:n-1
          M(i,i-1)=1;
          M(i,i)=-3*x(i)^2;
          end
      M(n,n-1)=1;
      M(n,n)=-1;
      end
%
%
%z = y = x = linspace (0, 1, 20);
%msh = bim3c_mesh_properties (msh3m_structured_mesh (x, y, z, 1, 1:6));
%x = msh.p(1, :)';
%y = msh.p(2, :)';
%z = msh.p(3, :)';
%
%nnodes = columns (msh.p);
%dnodes = bim3c_unknowns_on_faces (msh, [1, 3]);
%inodes = setdiff (1:nnodes, dnodes);
%
%A = bim3a_laplacian (msh, 1, 1);
%M = bim3a_reaction (msh, 1, 1);
%M  = A(inodes, inodes) + diag (diag (M(inodes, inodes)) .* (exp (u(inodes)) + exp (-u(inodes))));
%
%
%        end