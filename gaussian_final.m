clearvars;
tic;
K=3;
G1=[1 1 1];
G2=[1 0 1];
k=input('enter k: ');
u=generate_in(k);
original_msg=u;
u=[u 0 0];
rate=1/2;

%probability
SNRdb=0:0.5:8;
SNRlin=10.^(SNRdb/10);
p=qfunc(sqrt(2.*rate*SNRlin));

%encoding the message
op_en=encoder(k,G1,G2,u);

%output and new state tables
op_table=[0 3; 3 0; 2 1; 1 2];
ns_table=[0 2;0 2;1 3;1 3];

Nsim=20000;

for y=1:17
    Nerr=0;
    for x=1:1:Nsim
        %passing through bsc channel
        op_gaussian=gaussian_channel(op_en,rate,SNRlin(y),k);
        %generate trellis
        trellis=gen_trellis(op_gaussian,ns_table,op_table);
        %traceback
        op_gaussian_de=traceback(op_gaussian,trellis,ns_table,op_table);
        %calculate error
        err=cal_Nerr(k,u,op_gaussian_de);
        %update Net error
        Nerr=Nerr+err;
    end
    BER(y)=Nerr/(k*Nsim);
end
figure();
semilogy(SNRdb,BER,'o-','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
xlabel('SNR per Bit in dB'); 
ylabel('Probability of Bit Error'); 
grid on; 
toc;

%generate input
function [y] = generate_in(k)
    a=-ones(1,k);
    for i=1:1:k
        a(i)=randi(2)-1;
    end
    y=a;   
end

%gaussian channel
function [y]=gaussian_channel(op_en,rate,SNRlin,k)
    sigma = sqrt(1/(2*rate*SNRlin));
    for i = 1:2*(k+2)
        n = randn*sigma;
        op_en(i) = (2*op_en(i)-1) + n;
    end
    y=op_en;
end

%encoder
function [z]=encoder(k,G1,G2,u)
  op_en=zeros(0,0);
  R=zeros(1,3);
    for i=1:k+2
      y=u(i);
      R(3)=R(2);%moving digits in register forward
      R(2)=R(1);
      R(1)=y;%inserting 1 bit of input in the register
      B1=G1.*R;
      B2=G2.*R;
      op_en=[op_en,mod((B1(1)+B1(2)+B1(3)),2),mod((B2(1)+B2(2)+B2(3)),2)]; %updating the encoded array
    end
    z=op_en;
end

%calculate euclid distance
function [HD]=euclid_dist(A,B)
    HD=0;
    b=[];
    if(B==0)
        b=[-1,-1];
    elseif(B==1)
        b=[-1,1];
    elseif(B==2)
        b=[1,-1];
    elseif(B==3)
        b=[1,1];
    end
    HD=HD+(A(1)-b(1))*(A(1)-b(1))+(A(2)-b(2))*(A(2)-b(2));
end

function [y] = gen_trellis(ip,ns_table,op_table)
 
 trellis=-ones(4,(length(ip)/2)+1);
 trellis(1,1)=0;
 
 for i=1:length(ip)/2
   for j=1:4
     if(trellis(j,i)==-1) %if cell value is -1 then that cell is not included in any path
        trellis(j,i)=55;
        continue;
     end
     for k=1:2
       t = euclid_dist([ip(2*i-1),ip(2*i)],op_table(j,k));%comparison of pairwise input ip and output table
       if(trellis(ns_table(j,k)+1,i+1) == -1)
       trellis(ns_table(j,k)+1,i+1) = t+trellis(j,i);%insertion of hamming distance in the target cell
       else
       trellis(ns_table(j,k)+1,i+1) = min(t+trellis(j,i),trellis(ns_table(j,k)+1,i+1));%insertion of hamming distance in the target cell
     end     
     end
   end
 end
 y=trellis;
end

%ip is the input to decoder
%ns_table is the transition table from one state to another
%op_table is the output table for inputs from one state to another
function [y] = traceback(ip,trellis,ns_table,op_table)
    op_de=[];
    sz=(size(ip,2)/2)+1;
    last_col=[trellis(1,sz),trellis(2,sz),trellis(3,sz),trellis(4,sz)];
    last_col_min=min(last_col); %minimum of last column of the trellis
    state=0;
    for i=1:1:4
        if(last_col(i)==last_col_min)
            state=i; %state having minimum path metric
            break;
        end
    end
    for j=1:1:size(ip,2)/2
        state_previous=[];
        for x=1:1:2
            for y=1:1:4
                if(ns_table(y,x)==state-1)
                    state_previous=[state_previous,y-1]; %array containing the 2 possible previous states
                end
            end
        end
        input=[];
        for x=1:1:2
            for y=1:1:4
                if(ns_table(y,x)==state-1)
                    input=[input,x-1]; 
                    %array containing the inputs to go from previous state to the current state
                end
            end
        end
        sz=sz-1;
            str=[ip(2*sz-1),ip(2*sz)];
            str1=op_table(state_previous(1)+1,input(1)+1);
            str2=op_table(state_previous(2)+1,input(2)+1);
            an1 = hamm_dist(str,str1);
            an2 = hamm_dist(str,str2);
            if(an1+trellis(state_previous(1)+1,sz)>an2+trellis(state_previous(2)+1,sz))
                a=2;
            else
                a=1;
            end
        if(a==1)
            state=state_previous(1)+1;
        else
            state=state_previous(2)+1;
        end
        op_de=[input(a),op_de];
    end
    y=op_de;
    return;
end

function [y] = cal_Nerr(k,encoder_input,decoder_output)
    nerr=0;
    for i=1:k
        if(decoder_output(i) ~= encoder_input(i))%comparing decoded output and original message bits
            nerr=nerr+1;
        end
    end
    y=nerr;
    return;
end