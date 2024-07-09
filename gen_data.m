uii_range = linspace(3e-3,5e-3,10);
uij_range = linspace(0.5e-3,1.5e-3,10);
t_c_range = linspace(-1e-4,-1e-4,1);
t_so_range = linspace(1e-4,1e-4,1);
t_vo_range = linspace(0,0,1);
gs_range = linspace(1,10,10);
gv_range = linspace(10,40,10);
b_par_range = linspace(80e-6,80e-6,1);
b_perp_range = linspace(1350e-3,1350e-3,1);

k=1;

for uii = uii_range
    for uij = uij_range
        for t_c = t_c_range
            for t_so = t_so_range
                for t_vo = t_vo_range
                    for gs = gs_range
                        for gv = gv_range
                            for b_par = b_par_range
                                for b_perp = b_perp_range
                                    params = [uii;uij;t_c;t_so;t_vo;gs;gv;b_par;b_perp];
                                    k
                                    writematrix(params, "Indata/"+sprintf('%05d', k)+".csv"); 
                                    k=k+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end