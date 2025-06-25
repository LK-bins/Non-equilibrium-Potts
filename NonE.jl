# ndt (缺陷密度)
function latticesq(ll)
    nn = ll * ll  # 总格点数
    nbor = Array{Int}(undef, 8, nn)  # 第1-4列为最近邻，第5-8列为次近邻
    for s0 = 1:nn
        x0 = mod(s0 - 1, ll)
        y0 = div(s0 - 1, ll)
        
        # 最近邻坐标（右、上、左、下）
        x_right = mod(x0 + 1, ll)
        x_left  = mod(x0 - 1 + ll, ll)
        y_up    = mod(y0 + 1, ll)
        y_down  = mod(y0 - 1 + ll, ll)
        
        # 次近邻坐标（右上、右下、左上、左下）
        x_right_up   = mod(x0 + 1, ll)
        y_right_up   = mod(y0 + 1, ll)
        x_right_down = mod(x0 + 1, ll)
        y_right_down = mod(y0 - 1 + ll, ll)
        x_left_up    = mod(x0 - 1 + ll, ll)
        y_left_up    = mod(y0 + 1, ll)
        x_left_down  = mod(x0 - 1 + ll, ll)
        y_left_down  = mod(y0 - 1 + ll, ll)
        
        # 最近邻索引
        nbor[1, s0] = 1 + x_right + y0 * ll          # 右
        nbor[2, s0] = 1 + x0 + y_up * ll              # 上
        nbor[3, s0] = 1 + x_left + y0 * ll           # 左
        nbor[4, s0] = 1 + x0 + y_down * ll           # 下
        
        # 次近邻索引
        nbor[5, s0] = 1 + x_right_up + y_right_up * ll    # 右上
        nbor[6, s0] = 1 + x_right_down + y_right_down * ll # 右下
        nbor[7, s0] = 1 + x_left_up + y_left_up * ll      # 左上
        nbor[8, s0] = 1 + x_left_down + y_left_down * ll  # 左下
    end
    return nbor
end

 
# Constructs the initial random spin configuration
  function initspin(nn::Int)
     spin=Array{Int32}(undef,nn)
     for i=1:nn
        spin[i]=rand(1:1:3)
     end
     return spin
  end



function try_value(n, m)
tryvalue = Array{Int32}(undef, n, m)
tryvalue[1, 1] = 2; tryvalue[1, 2] = 3;
tryvalue[2, 1] = 1; tryvalue[2, 2] = 3;
tryvalue[3, 1] = 1; tryvalue[3, 2] = 2;
return tryvalue
end

# MC
function metropolis(ll::Int, nbor, spin, J1,J2, parT)
    nn::Int = ll * ll
    lx = ll
    ly = ll
    qn = 3      # q 值     # 邻居数
    t_value = try_value(3, 2)  # 确保 try_value 已定义

    for ix in 1:lx
        for iy in 1:ly
            ns = (ix - 1) * lx + iy
            rn = rand()
            value = min(floor(Int, rn * (qn - 1)) + 1, qn - 1)
            tr_value = t_value[spin[ns], value]

            E =  0
            E1 = 0
            for k in 1:4
                if spin[ns] == spin[nbor[k, ns]]
                    E += 1
                end
                if tr_value == spin[nbor[k, ns]]
                    E1 += 1
                end
            end
            
            E_nearest=0
            E1_nearest=0
            for k in 5:8
            if spin[ns] == spin[nbor[k, ns]]
                    E_nearest += 1
                end
                if tr_value == spin[nbor[k, ns]]
                    E1_nearest += 1
                end
            end
            
E_old=-J1*(E)-J2*(E_nearest)
E_new=-J1*(E1)-J2*(E1_nearest)
sum_delta_E=(E_new-E_old)    
           
            if parT>0
            prob=exp(-sum_delta_E/parT)
            elseif sum_delta_E<=0
                prob=1
            else
                prob=0
            end
            if rand() < prob  
                spin[ns] = tr_value
              end
        end
    end
    return spin
end


function energy(ll::Int, nbor, spin, J1, J2)
  nx=ll;
  nn=ll*ll;
 current_E = 0.0
for i = 1:nx*nx
                 sum_spin_nbor = 0
                 sum_spin_next_nbor = 0                        
                    for knbor = [1,2]                         #最近邻
                             id_nbor = nbor[knbor,i]
                                 if spin[id_nbor] == spin[i]
                                      sum_spin_nbor += 1
                                 end
                    end                
                    for knbor = [5,6]                          #次近邻
                             id_nbor = nbor[knbor,i]
                                 if spin[id_nbor] == spin[i]
                                     sum_spin_next_nbor += 1
                                 end
                    end                  
             current_E += -J1 * sum_spin_nbor - J2 * sum_spin_next_nbor    # 能量总和
                end
    return current_E/nn
end
  


@time begin
using Statistics
using Base.Threads  
nx=256
J1=-2;J2=1;
nbor=latticesq(nx)
Tc=2.035
m=11
#bins=1000
istep=1       
iistep=1000   # 样本序列
    for m in 5:5
        Tq=[10 20 40 80 160 320 640 1280 2560 5120 10240][m]
        parT_values = range(1.6,8.0, length=2*Tq)
     #    kinks = zeros(Float64,n_sample, 2*Tq)
         Energy= zeros(Float64,iistep,2*Tq)
         mss = zeros(Float64,iistep,2*Tq)
         ms2 = zeros(Float64,iistep,2*Tq)
         ms4 = zeros(Float64,iistep,2*Tq) 
     #    kink_avg=zeros(Float64,m,2*Tq)
         Energy_avg = zeros(Float64,2*Tq)
         ms_avg=zeros(Float64,2*Tq)
         ms2_avg=zeros(Float64,2*Tq)
         ms4_avg=zeros(Float64,2*Tq)
         #local_Energy = [zeros(Float64, 2*Tq) for _ in 1:nthreads()] 
         
      #   lock_obj = ReentrantLock()
        
    Base.Threads.@threads for sample_step = 1:iistep
            
         spin = initspin(nx*nx) 
         local_Energy =  zeros(Float64,2*Tq)   
           for t = 1:2*Tq
                 parT=9.6-parT_values[t] 
                        for n = 1:istep
                                spin = metropolis(nx, nbor, spin, J1,J2, parT)
                        end
                                current_E = energy(nx, nbor, spin, J1, J2)  
                                local_Energy[t]=current_E
           end  
                   Energy[sample_step,:]=local_Energy
     end  #  iistep
    
        for t=1:2*Tq
    Energy_avg[t]= mean(Energy[1:iistep, t])
        end
print("LK  Energy_avg Tq=:",Tq,'\n')
print(Energy_avg[:],'\n')         
end # m   
#savefig("D:\\2-Ising model\\Nd\\sample=100 nx=1024.png") 
end
