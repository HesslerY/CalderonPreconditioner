% D = d_mat_new(obj)
% Integration on testing functions: 1 point at cent
%
% Imput:
% object = RWG model (see object.m)
%
% Juan M. Rius, January 1997, modified Jan 1998 for struct object

function D = d_mat_new(obj)

Ne = size(obj.edges,2);
Nt = size(obj.trian,2);
D = sparse(Ne,Ne);

disp('computing D...');

v1 = obj.vertex(:,obj.topol(1,:));
v2 = obj.vertex(:,obj.topol(2,:));
v3 = obj.vertex(:,obj.topol(3,:));

r31 = v1-v3; 
r32 = v2-v3;
rs = r31+r32;

r(1,:,:)=r31;
r(2,:,:)=r32;
r(3,:,:)=zeros(size(r31));

l31 = sum(r31.^2,1);
l32 = sum(r32.^2,1);

fac1 = l31 + l32 + sum(r31.*r32,1);

for n=1:3,
   for m=1:3,
      
      en = obj.trian(n,:);
      em = obj.trian(m,:);
      s=sign(en).*sign(em);
      en=abs(en);
      em=abs(em);
      nm = en & em;
      
      tmp1 = 6*sum( squeeze(r(n,:,:)) .* squeeze(r(m,:,:)) ,1);
      tmp2 = 2*sum( squeeze(r(n,:,:)).* rs ,1);
      tmp3 = 2*sum( rs.* squeeze(r(m,:,:)) ,1);
      
      tmp1 = s(nm)/12 .* (fac1(nm) + tmp1(nm) - tmp2(nm) - tmp3(nm)) .* obj.ln(en(nm)) .* obj.ln(em(nm));
      if not(isempty(tmp1)),
      	tmp = tmp1 ./ (2*obj.ds(nm));
      	D = D + sparse(en(nm),em(nm),tmp,Ne,Ne);
      end
   end
end

