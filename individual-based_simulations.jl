using Random
using StatsBase
using Distributions
using DelimitedFiles

#parameter values
b::Float64=1.
u::Float64=0.45
cd::Float64=0.15
n::Int = 4
delta::Float64=3.
cost::Float64=0.3
n_traits::Int=n+2
n_deme::Int=4000
pmut::Float64=0.01
sigmut::Float64=0.01
gen::Int=150000
dirname="results"
sfreq::Int=1000


#initialization
pop_all=fill(0.0,(n+2,n,n_deme))
text0=string(0)
path=string(dirname,"/n_",n,"_u_",u,"_b_",b,"_cost_",cost,"_delta_",delta,"/cd_",cd )
mkpath(string(path,"/sampling")) #creates a folder for the output

for id_trait in 1:n_traits #allows user to enter initial values for each trait
    if id_trait<n_traits
        println(string("Enter initial value for dispersal with ",id_trait-1," cooperators"))
        ini_z=readline()
        ini_z = parse(Float64, ini_z)
    else
        println("Enter initial value for help ")
        ini_z=readline()
        ini_z = parse(Float64, ini_z)
    end
    pop_all[id_trait,:,:].=ini_z
    global text0=string(text0," ", ini_z, " ",0)
end

open(string(path,"/variables.txt"), "w") do io #creates file to keep track of average and variance in genetic traits
    println(io, text0)
end

namepara=string("b u help_cost disp_cost pmut sigmut nindiv ndeme ngen sampfreq") #creates file to keep track of the parameter used
allpara=string(b," ",u," ",cost," ",cd," ",pmut," ",sigmut," ",n," ",n_deme," ",gen," ",sfreq)
open(string(path,"/parameters.txt"), "w") do io
    println(io, namepara)
    println(io, allpara)
end



function allhelper(indiv_pop::Matrix{Float64},n::Int,n_deme::Int)
    #function that determines which individual in the population help or defects
    help::Matrix{Float64}=rand(Float64, (n, n_deme))
    help=help.<indiv_pop
    return help
end

 function investment(b::Float64,u::Float64,n_helper::Vector{Int},n::Int)
     #function that computes the benefits from helping in each patch (depending on the number of helpers)
    val::Vector{Float64}=1 .+delta*b.*(n_helper./n).^u
    return val
end

function pay_cost(pop_helping::Matrix{Float64},cost::Float64,delta::Float64)
    #cfunction that omputes the cost paid by each individual 
    val::Matrix{Float64}=delta.*cost.*pop_helping
    return val
end

n_helper=fill(0,n_deme)
disp_realized=fill(0.0,(n,n_deme))
popnext=fill(0.0,(n_traits,n,n_deme))
p::Matrix{Float64}=fill(0.0,(n,n_deme))

@time for it in 1:gen
    pop_helping=allhelper(pop_all[n_traits,:,:],n,n_deme) #determine which inividual helps
    for i in 1:n_deme
        n_helper[i]=Int(sum(pop_helping[:,i])) #assign to each patch its number of helper
        disp_realized[:,i].=pop_all[n_helper[i]+1,:,i] #assign to each individual its dispersal probability
    end
    fec::Matrix{Float64}=investment(b,u,n_helper,n)'.-pay_cost(pop_helping,cost,delta) #compute fecundity
    
    pop_vec=reshape(pop_all,(n_traits,n_deme*n))
    mutations=fill(0.0,(n_traits,n_deme*n))

    p[:,1]=fec[:,1].*(1 .-disp_realized[:,1]) 
    p[:,2:n_deme]=fec[:,2:n_deme].*(1-cd).*(disp_realized[:,2:n_deme]/(n_deme-1)) #compute the probability for each individual to have successful offspring in the first patch
    
    idfec=sample(1:(n*n_deme), Weights(vec(p)), n) #sample the offspring generation from the first patch
    popnext[:,:,1] = pop_vec[:,idfec]
    for i in 2:n_deme 
      p[:,i-1]=fec[:,i-1].*(1-cd).*(disp_realized[:,i-1]/(n_deme-1)) #update the probability of having successful offspring in the next patch 
      p[:,i]=fec[:,i].*(1 .-disp_realized[:,i])
      local idfec=sample(1:(n*n_deme), Weights(vec(p)), n) #sample the offspring generation in that patch
      popnext[:,:,i]=pop_vec[:,idfec]
    end
    idmut=sample(1:(n_deme*n), rand(Binomial(n*n_deme,pmut)),replace=false) #sample individuals that mutate
    mutations[:,idmut].=mutations[:,idmut].+rand(Normal(0,sigmut), (n_traits,length(idmut))) #sample their phenotypic perturbation
    global pop_all=popnext+reshape(mutations,(n_traits,n,n_deme)) #update their traits
    pop_all[pop_all.<0].=0 
    pop_all[pop_all.>1].=1 #truncate the traits to keep them between 0 and 1
    text=string(it)

        for id_trait in 1:n_traits #compute and record the average and variance in all traits
            mean_trait::Float64=sum(pop_all[id_trait,:,:])/(n_deme*n)
            var_trait::Float64=sum(pop_all[id_trait,:,:].*pop_all[id_trait,:,:])/(n_deme*n)-mean_trait^2
            text=string(text," ", mean_trait," ",var_trait)
        end
        open(string(path,"/variables.txt"), "a") do io
            println(io, text)
        end

    if it%sfreq==0 #every sfreq generation, sample the entire population
        println(it)
        for id_trait in 1:n_traits
            writedlm(string(path,"/sampling/x_",id_trait,"_generation_",it,".txt"), vec(pop_all[id_trait,:,:]))
        end
    end
end
  


 
