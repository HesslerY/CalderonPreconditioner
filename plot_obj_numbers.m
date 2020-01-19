% plot_obj_numbers(obj)

function plot_obj_numbers(obj)

Ne = length(obj.ln);
Nt = size(obj.un,2);
Nv = size(obj.vertex,2);
L = sqrt(sum(obj.ds)/6);  % Order of magnitude of object size 


% Plot triangle numbers
stt = cellstr(string(1:Nt));
rstt = obj.cent + 0.02*L*obj.un;
text(rstt(1,:), rstt(2,:), rstt(3,:), stt)

% Plot edge numbers
Tp = obj.edges(1,:);
Tm = obj.edges(2,:);
Vp = obj.topol(:,Tp);
Vm = obj.topol(:,Tm);
ve = reshape(Vp(Vp~=repmat(obj.edges(3,:),3,1)),[2 Ne]);
vcent = (obj.vertex(:,ve(1,:)) + obj.vertex(:,ve(2,:)))/2;
rste = vcent + 0.02*L* (obj.un(:,Tp) + obj.un(:,Tm))/2;
ste = cellstr(string(1:Ne));
text(rste(1,:), rste(2,:), rste(3,:), ste, 'color', 'red')

% Plot vertex numbers
stv = cellstr(string(1:Nv));
text(obj.vertex(1,:), obj.vertex(2,:), obj.vertex(3,:), stv, 'color', 'blue')
