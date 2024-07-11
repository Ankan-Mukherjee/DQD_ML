clear all; close all
warning off verbose
path_in = "Indata/";
infiles = dir(path_in+"/*.csv");
% Use all cores available
if isempty(gcp('nocreate'))
    parpool(feature('numcores'));
end
path1 = "Figures/" ;
path2 = "Outputdata/" ;
% rmdir Figures s;
% mkdir(path1);

%% Input Files
% var_idx_list = 1:size(infiles, 1);
var_idx_list=1:3;

%% Read variables from custom files
for var_idx = var_idx_list
vars = readmatrix(path_in + infiles(var_idx).name);
uii = vars(1); uij = vars(2); t_c = vars(3); gb_s = vars(4); gb_v = vars(5);

Batch_num = var_idx;
figname = "DQD"+num2str(Batch_num,'%05.f');
%% Do the computation
Vsd = -1500e-6;
n_ep = 100;
dVsd = Vsd*1e-2;
Vol= [Vsd-dVsd; Vsd + dVsd] ;
hbar=1.055e-34;q=1.609e-19;I0=q^2/hbar;muB=5.788e-5;muf = 0.0055;
kBT=5*0.0000156;

photo_relax = true;
section_cut = false;
row_sidx = idivide(int16(n_ep), 2); col_sidx = idivide(int16(n_ep), 1.7);

Hamiltonian = 0;

alpha21 = 0.5;alpha12 = 0.5;

eps1_max=5e-3;eps1_min = 1e-3;  % Combine both Forward and Rev Bias
eps2_max=-4e-3;eps2_min = -8.5e-3;

% Charging energy matrix Uij --> effect of electron addition in ith level on jth level
U = ones(4,4) * uij ;
U(1,1) = uii; U(2,2) = uii; U(3,3) = uii; U(4,4) = uii;
tt=t_c;  

gamma1=50e-6;                               % Lead coupling for lead 1
gamma2=50e-6;                               % Lead coupling for lead 2
factor = 0.03;                              % Photon relaxation magntude control

zp=gb_s * muB;                          % Zeeman splitting
vp=gb_v * muB;                         % Valley splitting

epsilon1_matrix = linspace(eps1_min, eps1_max, n_ep);
epsilon2_matrix = linspace(eps2_min, eps2_max, n_ep);
energy_matrix = zeros(n_ep,n_ep);
det_matrix = zeros(n_ep,n_ep);

numIterations = length(epsilon1_matrix);
% ppm = ParforProgressbar(numIterations);

I=zeros(n_ep, n_ep, 2);
G = zeros(n_ep, n_ep,1);
N_data = zeros(256, n_ep, n_ep, 2);
N1_data = zeros(256, n_ep, n_ep, 2);
N2_data = zeros(256, n_ep, n_ep, 2);


idx_mat = [2^7;2^6;2^5;2^4;2^3;2^2;2^1;2^0];
c1=zeros(256);c2=zeros(256);c3=zeros(256);c4=zeros(256);
c5=zeros(256);c6=zeros(256);c7=zeros(256);c8=zeros(256);


%   annihilation operators c1 ------> c8

for k1=1:2
    for k2=1:2
        for k3=1:2
            for k4=1:2
                for k5=1:2
                    for k6=1:2
                        for k7=1:2
                            kp = [1 k7-1 k6-1 k5-1 k4-1 k3-1 k2-1 k1-1] * idx_mat;
                            k =  [0 k7-1 k6-1 k5-1 k4-1 k3-1 k2-1 k1-1] * idx_mat;
                            c8(k+1,kp+1)=(-1)^(k7-1 + k6-1 + k5-1 + k4-1 + k3-1 + k2-1 + k1-1);
                        
                            kp = [k7-1 1 k6-1 k5-1 k4-1 k3-1 k2-1 k1-1] * idx_mat;
                            k =  [k7-1 0 k6-1 k5-1 k4-1 k3-1 k2-1 k1-1] * idx_mat;
                            c7(k+1,kp+1)=(-1)^(k6-1 + k5-1 + k4-1 + k3-1 + k2-1 + k1-1);
                        
                            kp = [k7-1 k6-1 1 k5-1 k4-1 k3-1 k2-1 k1-1] * idx_mat;
                            k =  [k7-1 k6-1 0 k5-1 k4-1 k3-1 k2-1 k1-1] * idx_mat;
                            c6(k+1,kp+1)=(-1)^(k5-1 + k4-1 + k3-1 + k2-1 + k1-1);

                            kp = [k7-1 k6-1 k5-1 1 k4-1 k3-1 k2-1 k1-1] * idx_mat;
                            k =  [k7-1 k6-1 k5-1 0 k4-1 k3-1 k2-1 k1-1] * idx_mat;
                            c5(k+1,kp+1)=(-1)^(k4-1 + k3-1 + k2-1 + k1-1);

                            kp = [k7-1 k6-1 k5-1 k4-1 1 k3-1 k2-1 k1-1] * idx_mat;
                            k =  [k7-1 k6-1 k5-1 k4-1 0 k3-1 k2-1 k1-1] * idx_mat;
                            c4(k+1,kp+1)=(-1)^(k3-1 + k2-1 + k1-1);

                            kp = [k7-1 k6-1 k5-1 k4-1 k3-1 1 k2-1 k1-1] * idx_mat;
                            k =  [k7-1 k6-1 k5-1 k4-1 k3-1 0 k2-1 k1-1] * idx_mat;
                            c3(k+1,kp+1)=(-1)^(k2-1 + k1-1);

                            kp = [k7-1 k6-1 k5-1 k4-1 k3-1 k2-1 1 k1-1] * idx_mat;
                            k =  [k7-1 k6-1 k5-1 k4-1 k3-1 k2-1 0 k1-1] * idx_mat;
                            c2(k+1,kp+1)=(-1)^(k1-1);

                            kp = [k7-1 k6-1 k5-1 k4-1 k3-1 k2-1 k1-1 1] * idx_mat;
                            k =  [k7-1 k6-1 k5-1 k4-1 k3-1 k2-1 k1-1 0] * idx_mat;
                            c1(k+1,kp+1)=1;
                        end
                    end
                end
            end
        end
    end
end


cn8 = c1;cn7 = c2;cn6 = c3;cn5 = c4;cn4 = c5;cn3 = c6;cn2 = c7;cn1 = c8;
c1 = cn1;c2 = cn2;c3 = cn3;c4 = cn4;c5 = cn5;c6 = cn6;c7 = cn7;c8 = cn8;

d1 = c1';d2 = c2';d3 = c3';d4 = c4';d5 = c5';d6 = c6';d7 = c7';d8 = c8'; %creation operators

del = 0.5e-3;

parfor ep1_idx = 1:n_ep
    for ep2_idx = 1:n_ep  
        for vsd_idx = 1:2    

            es1=-(epsilon1_matrix(ep1_idx) + alpha21 *epsilon2_matrix(ep2_idx)) ;
            es2=-(epsilon2_matrix(ep2_idx) + alpha12 *epsilon1_matrix(ep1_idx) + del);
            det_matrix(ep1_idx, ep2_idx) = es1 - es2;
            energy_matrix(ep1_idx, ep2_idx) = (es1 + es2);

            vec=[zeros(1,255) 1]';beta=0;
            H=zeros(256);N=zeros(256);N1 = zeros(256); N2 = zeros(256);    
%             zp = 0.1 * vp;
            % CONVENTION: 1st index is for the row, 2nd index is for the column

            ep=[(es1+zp+vp), (es1-zp+vp),es2+zp+vp, es2-zp+vp,...               % +zp --> up spin; +vp --> |+k> state 
                (es1+zp-vp), (es1-zp-vp),es2+zp-vp, es2-zp-vp];                 % -zp --> down spin; -vp --> |-k> state
            
            
            i8 = []; i7 = []; i6 = []; i5 = []; i4 = []; i3 = []; i2 = []; i1 = []; i0 = [];

            % diagonal terms of Hamiltonian matrix, H and number operator, N
            for k1=1:2
                for k2=1:2
                    for k3=1:2
                        for k4=1:2
                            for k5=1:2
                                for k6=1:2
                                    for k7=1:2
                                        for k8=1:2
                                            k = [k8-1 k7-1 k6-1 k5-1 k4-1 k3-1 k2-1 k1-1] * idx_mat;
                                            Ud=diag(diag(U));
                                            nv = [k7+k8-2;k5+k6-2;k3+k4-2; k1+k2-2];
                                            H(k+1,k+1)= (ep * [k8-1; k7-1; k6-1; k5-1; k4-1; k3-1; k2-1; k1-1]) + ...
                                                        (U(4,4)*((k1-1)&(k2-1)))+(U(3,3)*((k3-1)&(k4-1)))+...
                                                        (U(2,2)*((k5-1)&(k6-1)))+(U(1,1)*((k7-1)&(k8-1)))+...
                                                        0.5*((nv')*U*nv-(nv')*Ud*nv); 
                                            n = [k8-1 k7-1 k6-1 k5-1 k4-1 k3-1 k2-1 k1-1] * ones(8,1);
                                            N(k+1,k+1)= n;
                                            N2(k+1,k+1) = nv(1) + nv(3); N1(k+1, k+1) = nv(2) + nv(4);
                                            if n==0;  i0=[k+1;i0]; end
                                            if n==1;  i1=[k+1;i1]; end
                                            if n==2;  i2=[k+1;i2]; end
                                            if n==3;  i3=[k+1;i3]; end
                                            if n==4;  i4=[k+1;i4]; end
                                            if n==5;  i5=[k+1;i5]; end
                                            if n==6;  i6=[k+1;i6]; end
                                            if n==7;  i7=[k+1;i7]; end
                                            if n==8;  i8=[k+1;i8]; end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end


            % off-diagonal terms of [H] due to overlap

            % 1 -> L +K up
            % 2 -> L +K down
            % 3 -> R +K up
            % 4 -> R +K down
            % 5 -> L -K up
            % 6 -> L -K down
            % 7 -> R -K up
            % 8 -> R -K down

            for k1 = 1:2
                for k2 = 1:2
                    for k3 = 1:2
                        for k4 = 1:2
                            for k5 = 1:2
                                for k6 = 1:2
                                    % 1 <--> 3
                                    kp=[0, k6-1, 1, k5-1, k4-1, k3-1, k2-1, k1-1]*idx_mat;
                                    k=[1, k6-1, 0, k5-1, k4-1, k3-1, k2-1, k1-1]*idx_mat;
                                    H(k+1,kp+1)=tt*((-1)^(k6-1 + k4-1 + k2-1));
                                    H(kp+1,k+1)=H(k+1,kp+1);
                                    % 2 <--> 4
                                    kp=[k6-1, 0, k5-1, 1, k4-1, k3-1, k2-1, k1-1]*idx_mat;
                                    k=[k6-1, 1, k5-1, 0, k4-1, k3-1, k2-1, k1-1]*idx_mat;
                                    H(k+1,kp+1)=tt*((-1)^(k5-1 + k3-1 + k1-1));
                                    H(kp+1,k+1)=H(k+1,kp+1);
                                    % 5 <--> 7
                                    kp=[k6-1, k5-1, k4-1, k3-1, 0, k2-1, 1, k1-1]*idx_mat;
                                    k=[k6-1, k5-1, k4-1, k3-1, 1, k2-1, 0, k1-1]*idx_mat;
                                    H(k+1,kp+1)=tt*((-1)^(k6-1 + k4-1 + k2-1));
                                    H(kp+1,k+1)=H(k+1,kp+1);
                                    % 6 <--> 8
                                    kp=[k6-1, k5-1, k4-1, k3-1, k2-1, 0, k1-1, 1]*idx_mat;
                                    k=[k6-1, k5-1, k4-1, k3-1, k2-1, 1, k1-1, 0]*idx_mat;
                                    H(k+1,kp+1)=tt*((-1)^(k5-1 + k3-1 + k1-1));
                                    H(kp+1,k+1)=H(k+1,kp+1);
                                end
                            end
                        end
                    end
                end
            end

            Hamiltonian = H;

            %Separating the parts of the Hamiltonian with different numbers of electrons

            H0=H([i0],[i0]);H1=H([i1],[i1]);H2=H([i2],[i2]);
            H3=H([i3],[i3]);H4=H([i4],[i4]);H5=H([i5],[i5]);
            H6=H([i6],[i6]);H7=H([i7],[i7]);H8=H([i8],[i8]);
            
            N01 = N1(i0,i0);N11 = N1(i1,i1);N21 = N1(i2,i2);
            N31= N1(i3,i3); N41 = N1(i4, i4);N51 = N1(i5, i5);
            N61 = N1(i6, i6);N71 = N1(i7, i7);N81 = N1(i8, i8);

            N02 = N2(i0,i0);N12 = N2(i1,i1);N22 = N2(i2,i2);
            N32= N2(i3,i3); N42 = N2(i4, i4);N52 = N2(i5, i5);
            N62 = N2(i6, i6);N72 = N2(i7, i7);N82 = N2(i8, i8);

            %Finding eigenvalues and eigenvectors for Hamiltonians with 0,1,2,3and 4 electrons
            [V0,D0]=eig(H0);D0=real(sum(D0))'; [a0,b0]=sort(D0);
            [V1,D1]=eig(H1);D1=real(sum(D1))'; [a1,b1]=sort(D1);
            [V2,D2]=eig(H2);D2=real(sum(D2))'; [a2,b2]=sort(D2);
            [V3,D3]=eig(H3);D3=real(sum(D3))'; [a3,b3]=sort(D3);
            [V4,D4]=eig(H4);D4=real(sum(D4))'; [a4,b4]=sort(D4);
            [V5,D5]=eig(H5);D5=real(sum(D5))'; [a5,b5]=sort(D5);
            [V6,D6]=eig(H6);D6=real(sum(D6))'; [a6,b6]=sort(D6);
            [V7,D7]=eig(H7);D7=real(sum(D7))'; [a7,b7]=sort(D7);
            [V8,D8]=eig(H8);D8=real(sum(D8))'; [a8,b8]=sort(D8);


            %need annihilation operator between ’g’ and ’h’ subspaces
            c11=c1([i0],[i1]);c12=c2([i0],[i1]);c13=c3([i0],[i1]);c14=c4([i0],[i1]);
            c15=c5([i0],[i1]);c16=c6([i0],[i1]);c17=c7([i0],[i1]);c18=c8([i0],[i1]);

            c21=c1([i1],[i2]);c22=c2([i1],[i2]);c23=c3([i1],[i2]);c24=c4([i1],[i2]);
            c25=c5([i1],[i2]);c26=c6([i1],[i2]);c27=c7([i1],[i2]);c28=c8([i1],[i2]);

            c31=c1([i2],[i3]);c32=c2([i2],[i3]);c33=c3([i2],[i3]);c34=c4([i2],[i3]);
            c35=c5([i2],[i3]);c36=c6([i2],[i3]);c37=c7([i2],[i3]);c38=c8([i2],[i3]);

            c41=c1([i3],[i4]);c42=c2([i3],[i4]);c43=c3([i3],[i4]);c44=c4([i3],[i4]);
            c45=c5([i3],[i4]);c46=c6([i3],[i4]);c47=c7([i3],[i4]);c48=c8([i3],[i4]);

            c51=c1([i4],[i5]);c52=c2([i4],[i5]);c53=c3([i4],[i5]);c54=c4([i4],[i5]);
            c55=c5([i4],[i5]);c56=c6([i4],[i5]);c57=c7([i4],[i5]);c58=c8([i4],[i5]);

            c61=c1([i5],[i6]);c62=c2([i5],[i6]);c63=c3([i5],[i6]);c64=c4([i5],[i6]);
            c65=c5([i5],[i6]);c66=c6([i5],[i6]);c67=c7([i5],[i6]);c68=c8([i5],[i6]);

            c71=c1([i6],[i7]);c72=c2([i6],[i7]);c73=c3([i6],[i7]);c74=c4([i6],[i7]);
            c75=c5([i6],[i7]);c76=c6([i6],[i7]);c77=c7([i6],[i7]);c78=c8([i6],[i7]);

            c81=c1([i7],[i8]);c82=c2([i7],[i8]);c83=c3([i7],[i8]);c84=c4([i7],[i8]);
            c85=c5([i7],[i8]);c86=c6([i7],[i8]);c87=c7([i7],[i8]);c88=c8([i7],[i8]);

            %Evaluate #  and gamma2 using Fermi s Golden Rule
            s0=size([i0],1);s1=size([i1],1);s2=size([i2],1);s3=size([i3],1);s4=size([i4],1);
            s5=size([i5],1);s6=size([i6],1);s7=size([i7],1);s8=size([i8],1);

            eh1 = zeros(s0, s1);ev0 = zeros(1, s0);
            ev1 = zeros(s1, s1);
            for ct1=1:s0
                for ct2=1:s1
                    eh1(ct1,ct2)=a1(ct2)-a0(ct1);
                    ev0(:,ct1)=V0(:,b0(ct1));
                    ev1(:,ct2)=V1(:,b1(ct2));
%                     if abs(min(ev1(:,ct2))) == max(abs(ev1(:,ct2)))  
%                         ev1(:,ct2) = -1 * ev1(:,ct2);
%                     end
                end
            end

            A11=[c11*ev1]'*ev0;
            B11=[c12*ev1]'*ev0;
            C11=[c15*ev1]'*ev0;
            D11=[c16*ev1]'*ev0;
            gam11=gamma1*((A11.*A11)+(B11.*B11)+(C11.*C11)+(D11.*D11));
            A21=[c13*ev1]'*ev0;
            B21=[c14*ev1]'*ev0;
            C21=[c17*ev1]'*ev0;
            D21=[c18*ev1]'*ev0;
            gam21=gamma2*((A21.*A21)+(B21.*B21) + (C21.*C21)+(D21.*D21));

            eh2 = zeros(s1, s2); ev2 = zeros(s2, s2);
            for ct1=1:s1
                for ct2=1:s2
                    eh2(ct1,ct2)=a2(ct2)-a1(ct1);
                    ev2(:,ct2)=V2(:,b2(ct2));
%                     if abs(min(ev2(:,ct2))) == max(abs(ev2(:,ct2)))  
%                         ev2(:,ct2) = -1 * ev2(:,ct2);
%                     end
                end
            end
            A12=[c21*ev2]'*ev1;
            B12=[c22*ev2]'*ev1;
            C12=[c25*ev2]'*ev1;
            D12=[c26*ev2]'*ev1;
            gam12=gamma1*((A12.*A12)+(B12.*B12)+(C12.*C12)+(D12.*D12));
            A22=[c23*ev2]'*ev1;
            B22=[c24*ev2]'*ev1;
            C22=[c27*ev2]'*ev1;
            D22=[c28*ev2]'*ev1;
            gam22=gamma2*((A22.*A22)+(B22.*B22)+(C22.*C22)+(D22.*D22));

            eh3 = zeros(s2, s3); ev3 = zeros(s3, s3);
            for ct1=1:s2
                for ct2=1:s3
                    eh3(ct1,ct2)=a3(ct2)-a2(ct1);
                    ev3(:,ct2)=V3(:,b3(ct2));
%                     if abs(min(ev3(:,ct2))) == max(abs(ev3(:,ct2)))  
%                         ev3(:,ct2) = -1 * ev3(:,ct2);
%                     end
                end
            end
            A13=[c31*ev3]'*ev2;
            B13=[c32*ev3]'*ev2;
            C13=[c35*ev3]'*ev2;
            D13=[c36*ev3]'*ev2;
            gam13=gamma1*((A13.*A13)+(B13.*B13)+(C13.*C13)+(D13.*D13));
            A23=[c33*ev3]'*ev2;
            B23=[c34*ev3]'*ev2;
            C23=[c37*ev3]'*ev2;
            D23=[c38*ev3]'*ev2;
            gam23=gamma2*((A23.*A23)+(B23.*B23)+(C23.*C23)+(D23.*D23));

            eh4 = zeros(s3, s4); ev4 = zeros(s4, s4);
            for ct1=1:s3
                for ct2=1:s4
                    eh4(ct1,ct2)=a4(ct2)-a3(ct1);
                    ev4(:,ct2)=V4(:,b4(ct2));
%                     if abs(min(ev4(:,ct2))) == max(abs(ev4(:,ct2)))  
%                         ev4(:,ct2) = -1 * ev4(:,ct2);
%                     end
                end
            end
            A14=[c41*ev4]'*ev3;
            B14=[c42*ev4]'*ev3;
            C14=[c45*ev4]'*ev3;
            D14=[c46*ev4]'*ev3;
            gam14=gamma1*((A14.*A14)+(B14.*B14)+(C14.*C14)+(D14.*D14));
            A24=[c43*ev4]'*ev3;
            B24=[c44*ev4]'*ev3;
            C24=[c47*ev4]'*ev3;
            D24=[c48*ev4]'*ev3;
            gam24=gamma2*((A24.*A24)+(B24.*B24)+(C24.*C24)+(D24.*D24));

            eh5 = zeros(s4, s5); ev5 = zeros(s5, s5);
            for ct1=1:s4
                for ct2=1:s5
                    eh5(ct1,ct2)=a5(ct2)-a4(ct1);
                    ev5(:,ct2)=V5(:,b5(ct2));
%                     if abs(min(ev5(:,ct2))) == max(abs(ev5(:,ct2)))  
%                         ev5(:,ct2) = -1 * ev5(:,ct2);
%                     end
                end
            end
            A15=[c51*ev5]'*ev4;
            B15=[c52*ev5]'*ev4;
            C15=[c55*ev5]'*ev4;
            D15=[c56*ev5]'*ev4;
            gam15=gamma1*((A15.*A15)+(B15.*B15)+(C15.*C15)+(D15.*D15));
            A25=[c53*ev5]'*ev4;
            B25=[c54*ev5]'*ev4;
            C25=[c57*ev5]'*ev4;
            D25=[c58*ev5]'*ev4;
            gam25=gamma2*((A25.*A25)+(B25.*B25)+(C25.*C25)+(D25.*D25));

            eh6 = zeros(s5, s6); ev6 = zeros(s6, s6);
            for ct1=1:s5
                for ct2=1:s6
                    eh6(ct1,ct2)=a6(ct2)-a5(ct1);
                    ev6(:,ct2)=V6(:,b6(ct2));
%                     if abs(min(ev6(:,ct2))) == max(abs(ev6(:,ct2)))  
%                         ev6(:,ct2) = -1 * ev6(:,ct2);
%                     end
                end
            end
            A16=[c61*ev6]'*ev5;
            B16=[c62*ev6]'*ev5;
            C16=[c65*ev6]'*ev5;
            D16=[c66*ev6]'*ev5;
            gam16=gamma1*((A16.*A16)+(B16.*B16)+(C16.*C16)+(D16.*D16));
            A26=[c63*ev6]'*ev5;
            B26=[c64*ev6]'*ev5;
            C26=[c67*ev6]'*ev5;
            D26=[c68*ev6]'*ev5;
            gam26=gamma2*((A26.*A26)+(B26.*B26)+(C26.*C26)+(D26.*D26));

            eh7 = zeros(s6, s7); ev7 = zeros(s7, s7);
            for ct1=1:s6
                for ct2=1:s7
                    eh7(ct1,ct2)=a7(ct2)-a6(ct1);
                    ev7(:,ct2)=V7(:,b7(ct2));
%                     if abs(min(ev7(:,ct2))) == max(abs(ev7(:,ct2)))  
%                         ev7(:,ct2) = -1 * ev7(:,ct2);
%                     end
                end
            end
            A17=[c71*ev7]'*ev6;
            B17=[c72*ev7]'*ev6;
            C17=[c75*ev7]'*ev6;
            D17=[c76*ev7]'*ev6;
            gam17=gamma1*((A17.*A17)+(B17.*B17)+(C17.*C17)+(D17.*D17));
            A27=[c73*ev7]'*ev6;
            B27=[c74*ev7]'*ev6;
            C27=[c77*ev7]'*ev6;
            D27=[c78*ev7]'*ev6;
            gam27=gamma2*((A27.*A27)+(B27.*B27)+(C27.*C27)+(D27.*D27));

            eh8 = zeros(s7, s8); ev8 = zeros(s8, s8);
            for ct1=1:s7
                for ct2=1:s8
                    eh8(ct1,ct2)=a8(ct2)-a7(ct1);
                    ev8(:,ct2)=V8(:,b8(ct2));
                end
            end
            A18=[c81*ev8]'*ev7;
            B18=[c82*ev8]'*ev7;
            C18=[c85*ev8]'*ev7;
            D18=[c86*ev8]'*ev7;
            gam18=gamma1*((A18.*A18)+(B18.*B18)+(C18.*C18)+(D18.*D18));
            A28=[c83*ev8]'*ev7;
            B28=[c84*ev8]'*ev7;
            C28=[c87*ev8]'*ev7;
            D28=[c88*ev8]'*ev7;
            gam28=gamma2*((A28.*A28)+(B28.*B28)+(C28.*C28)+(D28.*D28));

            
            % ADD Photon Relaxation Terms
            eph1 = zeros(s1); eph2 = zeros(s2); eph3 = zeros(s3); eph4 = zeros(s4);
            eph5 = zeros(s5); eph6 = zeros(s6); eph7 = zeros(s7);

            Rph = zeros(256);

            % N = 1 space relaxations
            for ct1=1:s1
                for ct2=1:s1
                    eph1(ct1,ct2)=a1(ct1)-a1(ct2);
                end
            end
            a_ph = factor * gamma1;
            eph1(eph1==0)=NaN;
            b_eph1 = sign(eph1) .* a_ph .* (1./(exp(eph1/kBT)-1));
            b_eph1(isnan(b_eph1)|isinf(b_eph1))=0;

            % N = 2 space relaxations
            for ct1=1:s2
                for ct2=1:s2
                    eph2(ct1,ct2)=a2(ct1)-a2(ct2);
                end
            end
            a_ph = factor * gamma1;
            eph2(eph2==0)=NaN;
            b_eph2 = sign(eph2) .* a_ph .* (1./(exp(eph2/kBT)-1));
            b_eph2(isnan(b_eph2)|isinf(b_eph2))=0;

            % N = 3 space relaxations
            for ct1=1:s3
                for ct2=1:s3
                    eph3(ct1,ct2)=a3(ct1)-a3(ct2);
                end
            end
            a_ph = factor * gamma1;
            eph3(eph3==0)=NaN;
            b_eph3 = sign(eph3) .* a_ph .* (1./(exp(eph3/kBT)-1));
            b_eph3(isnan(b_eph3)|isinf(b_eph3))=0;

            % N = 4 space relaxations
            for ct1=1:s4
                for ct2=1:s4
                    eph4(ct1,ct2)=a4(ct1)-a4(ct2);
                end
            end

            a_ph = factor * gamma1;
            eph4(eph4==0)=NaN;
            b_eph4 = sign(eph4) .* a_ph .* (1./(exp(eph4/kBT)-1));
            b_eph4(isnan(b_eph4)|isinf(b_eph4))=0;

            % N = 5 space relaxations
            for ct1=1:s5
                for ct2=1:s5
                    eph5(ct1,ct2)=a5(ct1)-a5(ct2);
                end
            end

            a_ph = factor * gamma1;
            eph5(eph5==0)=NaN;
            b_eph5 = sign(eph5) .* a_ph .* (1./(exp(eph5/kBT)-1));
            b_eph5(isnan(b_eph5)|isinf(b_eph5))=0;

            % N = 6 space relaxations
            for ct1=1:s6
                for ct2=1:s6
                    eph6(ct1,ct2)=a6(ct1)-a6(ct2);
                end
            end

            a_ph = factor * gamma1;
            eph6(eph6==0)=NaN;
            b_eph6 = sign(eph6) .* a_ph .* (1./(exp(eph6/kBT)-1));
            b_eph6(isnan(b_eph6)|isinf(b_eph6))=0;

            % N = 7 space relaxations
            for ct1=1:s7
                for ct2=1:s7
                    eph7(ct1,ct2)=a7(ct1)-a7(ct2);
                end
            end

            a_ph = factor * gamma1;
            eph7(eph7==0)=NaN;
            b_eph7 = sign(eph7) .* a_ph .* (1./(exp(eph7/kBT)-1));
            b_eph7(isnan(b_eph7)|isinf(b_eph7))=0;

            Rph(s0+1:s0+s1,s0+1:s0+s1) = b_eph1;
            Rph(s0+s1+1:s0+s1+s2,s0+s1+1:s0+s1+s2) = b_eph2;
            Rph(s0+s1+s2+1:s0+s1+s2+s3,s0+s1+s2+1:s0+s1+s2+s3) = b_eph3;
            Rph(s0+s1+s2+s3+1:s0+s1+s2+s3+s4,s0+s1+s2+s3+1:s0+s1+s2+s3+s4) = b_eph4;
            Rph(s0+s1+s2+s3+s4+1:s0+s1+s2+s3+s4+s5,s0+s1+s2+s3+s4+1:s0+s1+s2+s3+s4+s5) = b_eph5;
            Rph(s0+s1+s2+s3+s4+s5+1:s0+s1+s2+s3+s4+s5+s6,s0+s1+s2+s3+s4+s5+1:s0+s1+s2+s3+s4+s5+s6) = b_eph6;
            Rph(s0+s1+s2+s3+s4+s5+s6+1:s0+s1+s2+s3+s4+s5+s6+s7,s0+s1+s2+s3+s4+s5+s6+1:s0+s1+s2+s3+s4+s5+s6+s7) = b_eph7;

            for de=1:256
                Rph(de,de)=-sum(Rph(:,de));
            end
            



            % Current and Rate Matrices
            vf=1;

            muL=muf(vf)+0.5*Vol(vsd_idx);
            muR=muf(vf)-0.5*Vol(vsd_idx);
            RL=zeros(256);
            R=zeros(256);

            % Fermi distribution parameters
            f11=1./(1+exp((eh1-muL)/kBT));
            f21=1./(1+exp((eh1-muR)/kBT));
            f12=1./(1+exp((eh2-muL)/kBT));
            f22=1./(1+exp((eh2-muR)/kBT));
            f13=1./(1+exp((eh3-muL)/kBT));
            f23=1./(1+exp((eh3-muR)/kBT));
            f14=1./(1+exp((eh4-muL)/kBT));
            f24=1./(1+exp((eh4-muR)/kBT));
            f15=1./(1+exp((eh5-muL)/kBT));
            f25=1./(1+exp((eh5-muR)/kBT));
            f16=1./(1+exp((eh6-muL)/kBT));
            f26=1./(1+exp((eh6-muR)/kBT));
            f17=1./(1+exp((eh7-muL)/kBT));
            f27=1./(1+exp((eh7-muR)/kBT));
            f18=1./(1+exp((eh8-muL)/kBT));
            f28=1./(1+exp((eh8-muR)/kBT));

            sti=1;stj=2;eni=s0;enj=1+s1;
            R(sti:eni,stj:enj)=gam11'.*(1-f11)+gam21'.*(1-f21);
            RL(sti:eni,stj:enj)=-gam11'.*(1-f11);
            R(stj:enj,sti:eni)=(gam11'.*f11+gam21'.*f21)';
            RL(stj:enj,sti:eni)=(gam11'.*f11)';

            sti=eni+1;stj=enj+1;eni=sti+s1-1;enj=stj+s2-1;
            R(sti:eni,stj:enj)=gam12'.*(1-f12)+gam22'.*(1-f22);
            RL(sti:eni,stj:enj)=-gam12'.*(1-f12);
            R(stj:enj,sti:eni)=(gam12'.*f12+gam22'.*f22)';
            RL(stj:enj,sti:eni)=(gam12'.*f12)';

            sti=eni+1;stj=enj+1;eni=sti+s2-1;enj=stj+s3-1;
            R(sti:eni,stj:enj)=gam13'.*(1-f13)+gam23'.*(1-f23);
            RL(sti:eni,stj:enj)=-gam13'.*(1-f13);
            R(stj:enj,sti:eni)=(gam13'.*f13+gam23'.*f23)';
            RL(stj:enj,sti:eni)=(gam13'.*f13)';

            sti=eni+1;stj=enj+1;eni=sti+s3-1;enj=stj+s4-1;
            R(sti:eni,stj:enj)=gam14'.*(1-f14)+gam24'.*(1-f24);
            RL(sti:eni,stj:enj)=-gam14'.*(1-f14);
            R(stj:enj,sti:eni)=(gam14'.*f14+gam24'.*f24)';
            RL(stj:enj,sti:eni)=(gam14'.*f14)';

            sti=eni+1;stj=enj+1;eni=sti+s4-1;enj=stj+s5-1;
            R(sti:eni,stj:enj)=gam15'.*(1-f15)+gam25'.*(1-f25);
            RL(sti:eni,stj:enj)=-gam15'.*(1-f15);
            R(stj:enj,sti:eni)=(gam15'.*f15+gam25'.*f25)';
            RL(stj:enj,sti:eni)=(gam15'.*f15)';

            sti=eni+1;stj=enj+1;eni=sti+s5-1;enj=stj+s6-1;
            R(sti:eni,stj:enj)=gam16'.*(1-f16)+gam26'.*(1-f26);
            RL(sti:eni,stj:enj)=-gam16'.*(1-f16);
            R(stj:enj,sti:eni)=(gam16'.*f16+gam26'.*f26)';
            RL(stj:enj,sti:eni)=(gam16'.*f16)';

            sti=eni+1;stj=enj+1;eni=sti+s6-1;enj=stj+s7-1;
            R(sti:eni,stj:enj)=gam17'.*(1-f17)+gam27'.*(1-f27);
            RL(sti:eni,stj:enj)=-gam17'.*(1-f17);
            R(stj:enj,sti:eni)=(gam17'.*f17+gam27'.*f27)';
            RL(stj:enj,sti:eni)=(gam17'.*f17)';

            sti=eni+1;stj=enj+1;eni=sti+s7-1;enj=stj+s8-1;
            R(sti:eni,stj:enj)=gam18'.*(1-f18)+gam28'.*(1-f28);
            RL(sti:eni,stj:enj)=-gam18'.*(1-f18);
            R(stj:enj,sti:eni)=(gam18'.*f18+gam28'.*f28)'; 
            RL(stj:enj,sti:eni)=(gam18'.*f18)';
            

            for de=1:255
                R(de,de)=-sum(R(:,de));
            end
            if (photo_relax == true) 
                R = R + Rph;
            end
            R(256,:)=1;
            if rank(R) < size(R,1)
                disp("Rank was:" + rank(R) + " idx1, idx2:" + ep1_idx +","...
                    + ep2_idx);
            end
            Pr=R\vec;
            Pr = double(Pr);
%             Pr = lsqminnorm(R,vec,1e-6);

%             tol = 1e-10;
%             maxit = size(R,1);Restart = size(R,1);
%             Pr = gmres(R,vec,Restart,tol,maxit);
% %             Pr = lsqr(R,vec,tol,1000);
%             Pr = Pr/sum(abs(Pr));
            I(ep1_idx, ep2_idx, vsd_idx)=sum(RL*Pr);          

            % Shift from eigenbasis to original tensor basis
            Pr_normal = zeros(256,1);
            Pr_normal(s0,1) =  Pr(s0,1);
            Pr_normal(s0+1:s0+s1,1) =  sum(Pr(s0+1:s0+s1,1)'.*ev1,2);
            Pr_normal(s0+s1+1:s0+s1+s2,1) =  sum(Pr(s0+s1+1:s0+s1+s2,1)'.*ev2,2);
            Pr_normal(s0+s1+s2+1:s0+s1+s2+s3,1) =  sum(Pr(s0+s1+s2+1:s0+s1+s2+s3,1)'.*ev3,2);
            Pr_normal(s0+s1+s2+s3+1:s0+s1+s2+s3+s4,1) =  sum(Pr(s0+s1+s2+s3+1:s0+s1+s2+s3+s4,1)'.*ev4,2);
            Pr_normal(s0+s1+s2+s3+s4+1:s0+s1+s2+s3+s4+s5,1) =  sum(Pr(s0+s1+s2+s3+s4+1:s0+s1+s2+s3+s4+s5,1)'.*ev5,2);
            Pr_normal(s0+s1+s2+s3+s4+s5+1:s0+s1+s2+s3+s4+s5+s6,1) =  sum(Pr(s0+s1+s2+s3+s4+s5+1:s0+s1+s2+s3+s4+s5+s6,1)'.*ev6,2);
            Pr_normal(s0+s1+s2+s3+s4+s5+s6+1:s0+s1+s2+s3+s4+s5+s6+s7,1) =  sum(Pr(s0+s1+s2+s3+s4+s5+s6+1:s0+s1+s2+s3+s4+s5+s6+s7,1)'.*ev7,2);
            Pr_normal(s0+s1+s2+s3+s4+s5+s6+s7+1:s0+s1+s2+s3+s4+s5+s6+s7+s8,1) =  sum(Pr(s0+s1+s2+s3+s4+s5+s6+s7+1:s0+s1+s2+s3+s4+s5+s6+s7+s8,1)'.*ev8,2);

            % Pr_normal = abs(Pr_normal);

            N_data(:,ep1_idx, ep2_idx, vsd_idx) = diag([0, ones(1,8), 2*ones(1,28), 3*ones(1,56),...
                4*ones(1,70), 5 * ones(1,56), 6*ones(1,28), 7*ones(1,8), 8]) * abs(Pr_normal);
            N1_data(:,ep1_idx, ep2_idx, vsd_idx)  =diag([diag(N01)', diag(N11)', diag(N21)', diag(N31)', diag(N41)',...
                                                            diag(N51)', diag(N61)', diag(N71)', diag(N81)']) * abs(Pr_normal);
            N2_data(:,ep1_idx, ep2_idx, vsd_idx) = diag([diag(N02)', diag(N12)', diag(N22)', diag(N32)', diag(N42)',...
                                                            diag(N52)', diag(N62)', diag(N72)', diag(N82)']) * abs(Pr_normal);
        end
    end
     pause(100/numIterations);
    % ppm.increment();
end
% delete(ppm);
G = (I(:, :, 2) - I(:, :, 1))/(2 * dVsd);
G = G /(q^2/(hbar * 2 * pi));

%% Plots for Current and conductance
I=I0*I; 
% max_I = max(max(I(:,:,2)));
% I(I < max_I*1e-2) = 0;
FigH = figure(108);
hold on
[x,y] = meshgrid(epsilon1_matrix,epsilon2_matrix);
% pcolor(y*1e3,x*1e3,log(abs(I(:, :, 2)) * 1e12));
pcolor(y*1e3,x*1e3,I(:, :, 2)' * 1e12);
h = colorbar;
% caxis([zMin, zMax]);
colormap jet
shading interp
set(gca,'FontSize',36)
xlabel('VG_2 (mV)','FontSize',36)
ylabel('VG_1 (mV)','FontSize',36)
title("DQD dot energy vs Current at V_{SD} = " + (Vsd)* 1e3 + " mV",'FontSize',40 )
% subtitle("I at t_c:"+abs(tt)+" U_{11}:"+U(1,1)+" B_{||}: "+B_par+" V_{SD}:"+Vsd+...
%     " B_{\perp}: "+B_perp + " tso: " + tso + " tvo: " + tvo , 'FontSize',20);
h.Label.String = "I_d (pA)";

if (section_cut == true)

    search_val = energy_matrix(row_sidx, col_sidx);
    [rows,cols] = find(abs(energy_matrix-search_val)< 1e-6);
    det_axis = diag(det_matrix(rows,cols));
    curr_det_axis = diag(I(rows, cols, 2) * 1e12);
    FigH = figure(72);
    plot(det_axis*(1e3),abs(curr_det_axis),'LineWidth',3);
    set(gca,'FontSize',48)
    xlabel('\delta (meV)','FontSize',48)
    ylabel('|I_D| (pA)','FontSize',48)
    yline(0,'--','y = 0','LineWidth',3)
    xline(-3,'--','x = -3','LineWidth',3)
    % subtitle("I at t_c:"+abs(tt)+" U_{11}:"+U(1,1)+" B_{||}: "+B_par+" V_{SD}:"+Vsd+...
    %     " B_{\perp}: "+B_perp + " tso: " + tso + " tvo: " + tvo + " \epsilon_{abs}: "+...
    %     search_val, 'FontSize',20);
    ylim padded;
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(FigH, path1+figname+'det.png','png');
    
end

FigH = figure(108);
if section_cut == true
% Draw a line along the detuning axis
    plot(epsilon2_matrix(cols)*1e3, epsilon1_matrix(rows)*1e3, '--', 'LineWidth', 4, 'Color', 'white');
end


set(FigH, 'Position', get(0, 'Screensize'));
saveas(FigH, path1+figname+'cmap.png','png');

end