clear all
close all
clc

syms vk akm gkm vm tk tm bkm bksh bkmsh a b c E F PG PGmin mpv

%%
fprintf('Derivadas Primeira e Segunda das Restrições de Igualdade e Desigualdade\n\n')

%% RESTRIÇÕES DE IGUALDADE %%

%% K é nó inicial

DPki = gkm*(vk/akm)^2-((vk*vm)/akm)*(gkm*cos(tk-tm)+bkm*sin(tk-tm));
DQki = -((bkm/(akm^2))+bkmsh)*vk^2+((vk*vm)/akm)*(bkm*cos(tk-tm)-gkm*sin(tk-tm));      

% Derivadas Primeira DeltaP
DPki_vk = diff(DPki,'vk');
DPki_vm = diff(DPki,'vm');
DPki_tk = diff(DPki,'tk');
DPki_tm = diff(DPki,'tm');
DPki_akm = diff(DPki,'akm');

% Derivadas Segunda DeltaPj
DPki_vkvk = diff(DPki_vk,'vk');
DPki_vkvm = diff(DPki_vk,'vm');
DPki_vmvm = diff(DPki_vm,'vm');

DPki_tkvk = diff(DPki_tk,'vk');
DPki_tmvk = diff(DPki_tm,'vk');
DPki_tkvm = diff(DPki_tk,'vm');
DPki_tmvm = diff(DPki_tm,'vm');

DPki_akmvk = diff(DPki_akm,'vk');
DPki_akmvm = diff(DPki_akm,'vm');

DPki_tktk = diff(DPki_tk,'tk');
DPki_tmtk = diff(DPki_tm,'tk');
DPki_tmtm = diff(DPki_tm,'tm');

DPki_akmtk = diff(DPki_akm,'tk');
DPki_akmtm = diff(DPki_akm,'tm');

DPki_akmakm = diff(DPki_akm,'akm');

% Derivadas Primeira DeltaQ
DQki_vk = diff(DQki,'vk');
DQki_vm = diff(DQki,'vm');
DQki_tk = diff(DQki,'tk');
DQki_tm = diff(DQki,'tm');
DQki_akm = diff(DQki,'akm');

% Derivadas Segunda DeltaQ
DQki_vkvk = diff(DQki_vk,'vk');
DQki_vkvm = diff(DQki_vk,'vm');
DQki_vmvm = diff(DQki_vm,'vm');

DQki_tkvk = diff(DQki_tk,'vk');
DQki_tmvk = diff(DQki_tm,'vk');
DQki_tkvm = diff(DQki_tk,'vm');
DQki_tmvm = diff(DQki_tm,'vm');

DQki_akmvk = diff(DQki_akm,'vk');
DQki_akmvm = diff(DQki_akm,'vm');

DQki_tktk = diff(DQki_tk,'tk');
DQki_tmtk = diff(DQki_tm,'tk');
DQki_tmtm = diff(DQki_tm,'tm');

DQki_akmtk = diff(DQki_akm,'tk');
DQki_akmtm = diff(DQki_akm,'tm');

DQki_akmakm = diff(DQki_akm,'akm');


%% K é nó final

DPkf = gkm*(vk^2)-((vk*vm)/akm)*(gkm*cos(tk-tm)+bkm*sin(tk-tm));
DQkf = -(bkm+bkmsh)*vk^2+((vk*vm)/akm)*(bkm*cos(tk-tm)-gkm*sin(tk-tm));

% Derivadas Primeira DeltaP
DPkf_vk = diff(DPkf,'vk');
DPkf_vm = diff(DPkf,'vm');
DPkf_tk = diff(DPkf,'tk');
DPkf_tm = diff(DPkf,'tm');
%DPkf_akm = diff(DPkf,'akm');
DPkf_akm = (vk*vm*(gkm*cos(tm - tk) + bkm*sin(tm - tk)))/akm^2;         %Utilizar essa equação para as derivadas de segunda ordem em relação aos taps/ângulos

% Derivadas Segunda DeltaP
DPkf_vkvk = diff(DPkf_vk,'vk');
DPkf_vkvm = diff(DPkf_vk,'vm');
DPkf_vmvm = diff(DPkf_vm,'vm');

DPkf_tkvk = diff(DPkf_tk,'vk');
DPkf_tmvk = diff(DPkf_tm,'vk');
DPkf_tkvm = diff(DPkf_tk,'vm');
DPkf_tmvm = diff(DPkf_tm,'vm');

DPkf_akmvk = diff(DPkf_akm,'vk');
DPkf_akmvm = diff(DPkf_akm,'vm');

DPkf_tktk = diff(DPkf_tk,'tk');
DPkf_tmtk = diff(DPkf_tm,'tk');
DPkf_tmtm = diff(DPkf_tm,'tm');

DPkf_akmtk = diff(DPkf_akm,'tk');
DPkf_akmtm = diff(DPkf_akm,'tm');

DPkf_akmakm = diff(DPkf_akm,'akm');

% Derivadas Primeira DeltaQ
DQkf_vk = diff(DQkf,'vk');
DQkf_vm = diff(DQkf,'vm');
DQkf_tk = diff(DQkf,'tk');
DQkf_tm = diff(DQkf,'tm');
%DQkf_akm = diff(DQkf,'akm');
DQkf_akm = -(vk*vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2;            %Utilizar essa equação para as derivadas de segunda ordem em relação aos taps/ângulos

% Derivadas Segunda DeltaQ
DQkf_vkvk = diff(DQkf_vk,'vk');
DQkf_vkvm = diff(DQkf_vk,'vm');
DQkf_vmvm = diff(DQkf_vm,'vm');

DQkf_tkvk = diff(DQkf_tk,'vk');
DQkf_tmvk = diff(DQkf_tm,'vk');
DQkf_tkvm = diff(DQkf_tk,'vm');
DQkf_tmvm = diff(DQkf_tm,'vm');

DQkf_akmvk = diff(DQkf_akm,'vk');
DQkf_akmvm = diff(DQkf_akm,'vm');

DQkf_tktk = diff(DQkf_tk,'tk');
DQkf_tmtk = diff(DQkf_tm,'tk');
DQkf_tmtm = diff(DQkf_tm,'tm');

DQkf_akmtk = diff(DQkf_akm,'tk');
DQkf_akmtm = diff(DQkf_akm,'tm');

DQkf_akmakm = diff(DQkf_akm,'akm');

%% RESTRIÇÕES DE DESIGUALDADE %%

%% K é nó inicial

Hki = -((bkm/akm^2)+bkmsh)*vk^2+((vk*vm)/akm)*(bkm*cos(tk-tm)-gkm*sin(tk-tm));

% Derivadas Primeira hi
Hki_vk = diff(Hki,'vk');
Hki_vm = diff(Hki,'vm');
Hki_tk = diff(Hki,'tk');
Hki_tm = diff(Hki,'tm');
Hki_akm = diff(Hki,'akm');

% Derivadas Segunda hi
Hki_vkvk = diff(Hki_vk,'vk');
Hki_vkvm = diff(Hki_vk,'vm');
Hki_vmvm = diff(Hki_vm,'vm');

Hki_tkvk = diff(Hki_tk,'vk');
Hki_tmvk = diff(Hki_tm,'vk');
Hki_tkvm = diff(Hki_tk,'vm');
Hki_tmvm = diff(Hki_tm,'vm');

Hki_akmvk = diff(Hki_akm,'vk');
Hki_akmvm = diff(Hki_akm,'vm');

Hki_tktk = diff(Hki_tk,'tk');
Hki_tktm = diff(Hki_tk,'tm');
Hki_tmtm = diff(Hki_tm,'tm');

Hki_akmtk = diff(Hki_akm,'tk');
Hki_akmtm = diff(Hki_akm,'tm');

Hki_akmakm = diff(Hki_akm,'akm');


%% K é no final

Hkf = -(bkm+bkmsh)*vk^2+((vk*vm)/akm)*(bkm*cos(tk-tm)-gkm*sin(tk-tm));

% Hkfneg = -Hkf;
% 
% Hkf_vk_neg = diff(Hkfneg,'vk');

% Derivadas Primeira hi
Hkf_vk = diff(Hkf,'vk');
Hkf_vm = diff(Hkf,'vm');
Hkf_tk = diff(Hkf,'tk');
Hkf_tm = diff(Hkf,'tm');
%Hkf_akm = diff(Hkf,'akm');
Hkf_akm = -(vk*vm*(bkm*cos(tm - tk) - gkm*sin(tm - tk)))/akm^2;             %Utilizar essa equação para as derivadas de segunda ordem

% Derivadas Segunda hi
Hkf_vkvk = diff(Hkf_vk,'vk');
Hkf_vkvm = diff(Hkf_vk,'vm');
Hkf_vmvm = diff(Hkf_vm,'vm');

Hkf_tkvk = diff(Hkf_tk,'vk');
Hkf_tmvk = diff(Hkf_tm,'vk');
Hkf_tkvm = diff(Hkf_tk,'vm');
Hkf_tmvm = diff(Hkf_tm,'vm');

Hkf_akmvk = diff(Hkf_akm,'vk');
Hkf_akmvm = diff(Hkf_akm,'vm');

Hkf_tktk = diff(Hkf_tk,'tk');
Hkf_tktm = diff(Hkf_tk,'tm');
Hkf_tmtm = diff(Hkf_tm,'tm');

Hkf_akmtk = diff(Hkf_akm,'tk');
Hkf_akmtm = diff(Hkf_akm,'tm');

Hkf_akmakm = diff(Hkf_akm,'akm');

%% FUNÇÃO OBJETIVO - MINIMIZAÇÃO DAS PERDAS

f = gkm*(((vk/akm)^2)+(vm^2)-2*(vk*vm/akm)*cos(tk-tm));

%Derivadas Primeira da Função Objetivo
f_vk = diff(f,'vk');
f_vm = diff(f,'vm');
f_tk = diff(f,'tk');
f_tm = diff(f,'tm');
f_akm = diff(f,'akm');

%Derivadas Segunda da Função Objetivo
f_vkvk = diff(f_vk,'vk');
f_vkvm = diff(f_vk,'vm');
f_vmvm = diff(f_vm,'vm');

f_tkvk = diff(f_tk,'vk');
f_tmvm = diff(f_tm,'vm');
f_tmvk = diff(f_tm,'vk');
f_tkvm = diff(f_tk,'vm');

f_akmvk = diff(f_akm,'vk');
f_akmvm = diff(f_akm,'vm');

f_vktk = diff(f_vk,'tk');
f_vmtm = diff(f_vm,'tm');
f_vmtk = diff(f_vm,'tk');
f_vktm = diff(f_vk,'tm');

f_tktk = diff(f_tk,'tk');
f_tktm = diff(f_tk,'tm');
f_tmtm = diff(f_tm,'tm');

f_akmtk = diff(f_akm,'tk');
f_akmtm = diff(f_akm,'tm');

f_akmakm = diff(f_akm,'akm');

%% FUNÇÃO OBJETIVO - CUSTO DE GERAÇÃO

fg = (a*PG^2+b*PG+c)/10000 + (E*mpv)/10000;

fg_PG = diff(fg,'PG');

fg_PGPG = diff(fg_PG,'PG');

fg_mpv = diff(fg,'mpv');

fg_mpv2 = diff(fg_mpv,'mpv');

%% RESTRIÇÃO DE DESIGUALDADE - PONTO DE VÁLVULA

%% Primera parte -h - mpv < 0

fpv1 = -(sin(F*(PGmin-PG)))-mpv;

fpv1_PG = diff(fpv1,'PG');

fpv1_PGPG = diff(fpv1_PG,'PG');

fpv1_m = diff(fpv1,'mpv');

fpv1_mm = diff(fpv1_m,'mpv');


%% Segunda parte h - mpv < 0

fpv2 = (sin(F*(PGmin-PG)))- mpv;

fpv2_PG = diff(fpv2,'PG');

fpv2_PGPG = diff(fpv2_PG,'PG');

fpv2_m = diff(fpv2,'mpv');

fpv2_mm = diff(fpv2_m,'mpv');





















