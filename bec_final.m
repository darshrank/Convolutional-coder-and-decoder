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

op_table=[0 3; 3 0; 2 1; 1 2];
ns_table=[0 2;0 2;1 3;1 3];

Nsim=20000;
for y=1:17
    Nerr=0;
    for x=1:1:Nsim
        %passing through bsc channel
        op_bec=bec_channel(op_en,k,p(y));
        ip_bec_de=[];
        for i=1:2:length(op_bec)
          if(op_bec(i)==0)
            if(op_bec(i+1)==0)
            ip_bec_de=[ip_bec_de,"00"];
            elseif(op_bec(i+1)==1)
            ip_bec_de=[ip_bec_de,"01"];
            elseif(op_bec(i+1)==5)
            ip_bec_de=[ip_bec_de,"05"];
            end
          elseif(op_bec(i)==1)
            if(op_bec(i+1)==0)
            ip_bec_de=[ip_bec_de,"10"];
            elseif(op_bec(i+1)==1)
            ip_bec_de=[ip_bec_de,"11"];
            elseif(op_bec(i+1)==5)
            ip_bec_de=[ip_bec_de,"15"];
            end
          else
            if(op_bec(i+1)==0)
            ip_bec_de=[ip_bec_de,"50"];
            elseif(op_bec(i+1)==1)
            ip_bec_de=[ip_bec_de,"51"];
            elseif(op_bec(i+1)==5)
            ip_bec_de=[ip_bec_de,"55"];
            end
          end
        end
        trellis=gen_trellis(ip_bec_de,ns_table,op_table);
        op_bsc_de=traceback(ip_bec_de,trellis,ns_table,op_table);

        err=cal_Nerr(k,u,op_bsc_de);
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

%bec channel
function [y]=bec_channel(op_en,k,p)
  for i=1:2*(k+2)
      b=rand(1);
      if(b<p)
          op_en(i)=5;
      end
  end
  y=op_en;
end

%encoder
function [z]=encoder(k,G1,G2,u)
  op_en=zeros(0,0);
  R=zeros(1,3);
    for i=1:k+2
      y=u(i);
      R(3)=R(2);
      R(2)=R(1);
      R(1)=y;
      B1=G1.*R;
      B2=G2.*R;
      op_en=[op_en,mod((B1(1)+B1(2)+B1(3)),2),mod((B2(1)+B2(2)+B2(3)),2)];
    end
    z=op_en;
end

%calculate hamming distance
function [HD]=hamm_dist(A,B)
HD=0;
if(B==0)
    B="00";
elseif(B==1)
    B="01";
elseif(B==2)
    B="10";
elseif(B==3)
    B="11";
end
for i=1:2
    if((A{1}(i)=='0' || A{1}(i)=='1') && (B{1}(i)=='0' || B{1}(i)=='1') )
      if(A{1}(i)~=B{1}(i))
       HD=HD+1;
      end
    else
      HD=HD+1;
    end
end
  
end

function [y] = gen_trellis(ip,ns_table,op_table)
 
 trellis=-ones(4,length(ip)+1);
 trellis(1,1)=0;
 
 for i=1:length(ip)
   for j=1:4
     if(trellis(j,i)==-1) %if cell value is -1 then that cell is not included in any path
        trellis(j,i)=55;
        continue;
     end
     for k=1:2
       t = hamm_dist(ip(i),op_table(j,k));%comparison of pairwise input ip and output table
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

    
function [y] = traceback(ip,trellis,ns_table,op_table)
    op_de=[];
    [a,c]=min(trellis,[],1);
    sz=size(ip,2)+1;
    a=a(sz);
    c=c(sz);
    for j=1:1:size(ip,2)
        bd=mod(find(ns_table==(c-1))-1,size(ns_table,1));
        bd=bd';
        ka=find(ns_table==(c-1))/size(ns_table,1);
        ka=ka';
        for i=1:1:size(ka,2)
            if(ka(i)<=1)
                ka(i)=0;
            else
                ka(i)=1;
            end
        end
        sz=sz-1;
            arr=ip(sz);
            arr1=op_table(bd(1)+1,ka(1)+1);
            arr2=op_table(bd(2)+1,ka(2)+1);
            an1 = hamm_dist(arr,arr1);
            an2 = hamm_dist(arr,arr2);
            if(an1+trellis(bd(1)+1,sz)>an2+trellis(bd(2)+1,sz))
                e=2;
            else
                e=1;
            end
        if(e==1)
            c=bd(1)+1;
        else
            c=bd(2)+1;
        end
        op_de=[ka(e),op_de];
    end
    y=op_de;
    return;
end

function [y] = cal_Nerr(k,encoder_input,decoder_output)
    nerr=0;
    for i=1:k
        if(decoder_output(i) ~= encoder_input(i))
            nerr=nerr+1;
        end
    end
    y=nerr;
    return;
end