using CairoMakie, DataFrames, LinearAlgebra, SparseArrays, Arpack, ProgressBars, IterTools, Combinatorics, DelimitedFiles, Tables, Random, Statistics
#commit?
println("\nstart")
"""N is the number of sites and M is the number of particles"""
N=10
M=5#Int64(floor(N/2))
parts=2
repuls=-10

pichart=collect(Float64, -pi:2*pi/N:(2*pi*(N-1)/N)-pi)
"""
potbonus adds a background potential energy
potweight is more important and determines the size of the random potential 1 is the same weight as the hopping term (plus or minus)
"""


""" Ton options are: "Hop","Hop+Rep","Hop+Pot","Hop+Pot+Rep" """
Ton="Hop+Pot"
alp=0
""" half can be 1 to get a projector that has two halfs or anything else for a staggered projector """
half=0
offset=0
doeik=0
bin=binomial(N,M)
start=0
ending=bin
"""
Aim lets you aim your half projector at the lowest potential in the anderson model
average lets you aim at the lowest average half
sizepar determines the size of the first half, the second half wil be N-sizepar
sitenump is the set of indivdual projectors (by default: first half has 0,1,2,3,....and M particles)
"""
normaly=0
aim=0
allop=1
average=1
sizepar=Int64(floor(N/2))
sitenump=Array(0:M)
if any(sitenump.<0)
    error("did you mean for sitenump to go below 0?
    If so comment out this error")
end

""" order determines to which power you want to take the partial overlaps """
minorder=1
maxorder=1
""" here a dictionary of youngtableaus is created """
name=zeros(1, M)
youngdict = Dict(1 => name[1:M])

print("whiles \n")
for i in 1:(bin-1)
    name[1]+=1

    while N-M+1 in name

        global place=findall(name.==N-M+1)
        name[place[end][2]+1]+=1
        name[1:place[end][2]+1].=name[place[end][2]+1]
    end
    youngdict[i+1]=name[1:M]

end
print("no more whiles ")
""" the inverse dictionairy of young tab """
inv_yd = Dict(values(youngdict) .=> keys(youngdict))
"""
translating the youngtableaus to locations of particles on the lattice
ttotal[k] gives the positions of the particles that correspond to to k-th entry in an eigenvector
"""
extra1=reverse(collect(0:M))
extra=extra1[1:M]
nametaart=copy(extra)
ttotal=Dict(1.0=>copy(nametaart))
for i in 2:bin
    ttotal[i]=copy(youngdict[i]+extra)
end
""" locfind does the opposite of ttotal keys are of the form [a,b,c,.......,z] of size M """
locfind=Dict(values(ttotal) .=> keys(ttotal))

"""
if halfpot is not 1 or 0 the system has a high but random potential of order 10 t

pot is the random potential per site
the lower 2 pots are some premade pots (one for a 6 site and one for a 12 site case) these can be used to compare.
It is also possible to comment out all pot before running the code again to preserve the last used pot

"""

#for bonusweight in [collect(0:0.1:1); collect(2:10)]
bonusweight=10
halfpot=1

if halfpot==1
    potter=[1:1:Int64(N/2);]
    bonusweight=bonusweight
    potweight=0.1
elseif halfpot==0
    potter=[1:2:N;]
    bonusweight=bonusweight
    potweight=0.1
elseif halfpot==2
    bonusweight=0
    potter=[1]
    potweight=10
elseif halfpot == 3
    bonusweight=0
    potter=[1]
    potweight=0.0001
end



potbonus=zeros(Float64,N)

potbonus[potter].+=bonusweight

currentseed=54321
Random.seed!(currentseed)

pot=potweight*(2*rand(Float64, N).-1)
pot.-=mean(pot)
pot.+=potbonus
pt=copy(pot)


""" Hp is the term in the hamiltonian that adds the anderson potential """

Hp=zeros(bin,bin)
for j in 1:bin
    valuetemp=0
    loctemp=extra+youngdict[j]
    for l in 1:M
        valuetemp+=pot[Int64(loctemp[l])]
    end
    Hp[j,j]=valuetemp
end

""" Htest is the hopping term and Hrep is the repulsion (or atraction if negative). Both are indepentitly created but in the same loop"""
jsneigh=[]
Htest=spzeros(bin,bin)
Hrep=spzeros(bin,bin)
if M > 1
    for j in 1:bin
        normy=0
        tablfor=copy(youngdict[j])
        if tablfor[1]!=N-M
            tablfor[1]+=1
            loc=inv_yd[tablfor]
            Htest[j,loc]=1
            normy+=1
        elseif tablfor[M]!=0
            tablfor[2:end].-=1
            tablfor[1]=0
            tablfor=circshift(tablfor,M-1)
            loc=inv_yd[tablfor]
            Htest[j,loc]=1
            normy+=1
        end
        if youngdict[j][1] == youngdict[j][2]
            Hrep[j,j]+=repuls
            push!(jsneigh,j)
        end
        for i in 2:M-1
            tabl=copy(youngdict[j])
            if tabl[i] < tabl[i-1]
                tabl[i]+=1
                loc=inv_yd[tabl]
                Htest[j,loc]=1
                normy+=1
            end
            if youngdict[j][i] == youngdict[j][i+1]
                Hrep[j,j]+=repuls
                push!(jsneigh,j)

            end
        end


        tablforpost=copy(youngdict[j])
        if tablforpost[M] < tablforpost[M-1]
            tablforpost[M]+=1
            loc=inv_yd[tablforpost]
            Htest[j,loc]=1
            normy+=1
        elseif tablforpost[M]==N-M
            tablforpost[1:M-1].-=1
            tablforpost[M]=0
            loc=inv_yd[tablforpost]
            Htest[j,loc]=1
            normy+=1
        end
        if youngdict[j][M] == 0.0 && youngdict[j][1]==Float64(N-M)
            Hrep[j,j]+=repuls
            push!(jsneigh,j)
            #print("did get here",j,Htest[j,j])

        end

        if normy!=0 && normaly==1
            Htest[j,:].*=1/sqrt(normy)
        end
    end
else
    for j in 1:bin


        tablfor=copy(youngdict[j])
        if tablfor[1]!=N-M
            tablfor[1]+=1
            loc=inv_yd[tablfor]
            Htest[j,loc]=1
        elseif tablfor[M]!=0
            tablfor[2:M].-=1
            tablfor[1]=0
            tablfor=circshift(tablfor,-1)
            loc=inv_yd[tablfor]
            Htest[j,loc]=1
        end
    end
end

Hhop=exp(im*alp/N)*Htest+exp(-im*alp/N)*transpose(Htest)
println("Hhop done")

if Ton=="Hop"
    Ham=Array(Hhop)
    potdis=""
elseif Ton=="Hop+Rep"
    Ham=Array(Hhop+Hrep)
    #push!(Qlist,Array(Hrep))
    potdis=""
elseif Ton=="Hop+Pot"
    Ham=Array(Hhop+Hp)
    #push!(Qlist,Array(Hp))
    potdis=string(" ",pot)
elseif Ton=="Hop+Pot+Rep"
    Ham=Array(Hhop+Hp+Hrep)
    # push!(Qlist,Array(Hp))
    # push!(Qlist,Array(Hrep))
    potdis=string(" ",pot)
else
    error("incorrect Ton")
end
println("ham done")

"""
Here the Qs are made. Which is essentialy a hop by a certain number of steps.
For now I make all the parts of the Qs explicitly even though Q_1 can be constructed from Q_1_right+h.c
and the Qs for current can be constructed from two parts of the other Qs as i*Q_1_right-i*Q_1_left for example.
"""

if Ton=="Hop"
    ccfor=zeros(ComplexF64,Int64(N*N/2),bin,bin)
    ccback=zeros(ComplexF64,Int64(N*N/2),bin,bin)
    Qlist=[]
    for q in 1:N/2
        """I do a loop twice. The extra one is for the h.c. This is not necessary but I wanted to be sure"""
        global Q_temp=zeros(bin,bin)
        for i in 1:bin
            for j in 1:N
                sign=1
                """ttotal is a list of positions of the particles on the lattice. ttotal[1] is [2,1] which represents (1100), ttotal[2] is [3,1] which represents (1010), etc."""
                tablfor=copy(ttotal[i])
                """if statement encodes the periodicity"""
                if j+q > N
                    #print("are")
                    qnew=q-N
                    #sign*=(-1)
                else
                    #print("we")
                    qnew=q
                end
                """
                If there is a particle at j and not at j+q replace the particle at j by a particle at j+q.
                locfind is the reverse dictionairy of ttotal locfind[[2,1]]=1 etc.
                """
                if (j+qnew) ∉ tablfor && (j) ∈ tablfor
                    tablfor[findfirst(tablfor .== j)] = (j+qnew)
                    global counter=0

                    while counter<(length(tablfor)-1)
                        global counter = 0
                        for k in 1:length(tablfor)-1
                            firstentry=tablfor[k]
                            secondentry=tablfor[k+1]
                            if firstentry<secondentry
                                tablfor[k]=secondentry
                                tablfor[k+1]=firstentry
                                sign*=-1
                            else
                                global counter +=1
                            end
                        end
                    end
                    ccfor[Int64((q-1)*N+j),Int64(locfind[tablfor]),Int64(locfind[copy(ttotal[i])])]=sign
                    loc=Int64(locfind[tablfor])
                    global Q_temp[loc,i]+=copy(sign)
                end
            end
            for j in 1:N
                sign=1
                """ttotal is a list of positions of the particles on the lattice. ttotal[1] is [2,1] which represents (1100), ttotal[2] is [3,1] which represents (1010), etc."""
                tablback=copy(ttotal[i])
                """if statement encodes the periodicity"""
                if j-q < 1
                    #print("doing")
                    qnew=q-N
                    #sign*=(-1)
                else
                    #print("this")
                    qnew=q
                end
                """
                If there is a particle at j and not at j-q replace the particle at j by a particle at j-q.
                """
                if (j-qnew) ∉ tablback && (j) ∈ tablback
                    tablback[findfirst(tablback .== j)] = (j-qnew)
                    global counter=0

                    while counter<(length(tablback)-1)
                        global counter = 0
                        for k in 1:length(tablback)-1
                            firstentry=tablback[k]
                            secondentry=tablback[k+1]
                            if firstentry<secondentry
                                tablback[k]=secondentry
                                tablback[k+1]=firstentry
                                sign*=-1
                            else
                                global counter +=1
                            end
                        end
                    end
                    ccback[Int64((q-1)*N+j),Int64(locfind[tablback]),Int64(locfind[copy(ttotal[i])])]=sign
                    loc=Int64(locfind[tablback])
                    global Q_temp[loc,i]+=copy(sign)
                end
            end
        end
        #"""transpose in this case is the same as the h.c."""
        Q_temp2=Q_temp#+transpose(Q_temp)
        push!(Qlist,copy(Q_temp2))
    end
    for q in 1:N/2-1
        global Q_temp=complex(zeros(bin,bin))
        for i in 1:bin
            for j in 1:N
                """This is where the i is introduced"""
                sign=im
                """ttotal is a list of positions of the particles on the lattice. In the N=4,M=2 case for example ttotal[1] is [2,1] which represents (1100), ttotal[2] is [3,1] which represents (1010), etc."""
                tablfor=copy(ttotal[i])
                """if statement encodes the periodicity"""
                if j+q > N
                    #print("are")
                    qnew=q-N
                    #sign*=(-1)
                else
                    #print("we")
                    qnew=q
                end
                """
                If there is a particle at j and not at j+q replace the particle at j by a particle at j+q.
                """
                if (j+qnew) ∉ tablfor && (j) ∈ tablfor
                    tablfor[findfirst(tablfor .== j)] = (j+qnew)
                    global counter=0

                    while counter<(length(tablfor)-1)
                        global counter = 0
                        for k in 1:length(tablfor)-1
                            firstentry=tablfor[k]
                            secondentry=tablfor[k+1]
                            if firstentry<secondentry
                                tablfor[k]=secondentry
                                tablfor[k+1]=firstentry
                                sign*=-1
                            else
                                global counter +=1
                            end
                        end
                    end
                    loc=Int64(locfind[tablfor])
                    global Q_temp[loc,i]+=copy(sign)
                end
            end
            for j in 1:N
                """ in the loop for the hermitian conjugate I use -i instead of i """
                sign=-im
                """ttotal is a list of positions of the particles on the lattice. ttotal[1] is [2,1] which represents (1100), ttotal[2] is [3,1] which represents (1010), etc."""
                tablback=copy(ttotal[i])
                """if statement encodes the periodicity"""
                if j-q < 1
                    #print("doing")
                    qnew=q-N
                    #sign*=(-1)
                else
                    #print("this")
                    qnew=q
                end
                """
                If there is a particle at j and not at j+q replace the particle at j by a particle at j+q.
                locfind is the reverse dictionairy of ttotal locfind[[2,1]]=1 etc.
                """
                if (j-qnew) ∉ tablback && (j) ∈ tablback
                    tablback[findfirst(tablback .== j)] = (j-qnew)
                    global counter=0

                    while counter<(length(tablback)-1)
                        global counter = 0
                        for k in 1:length(tablback)-1
                            firstentry=tablback[k]
                            secondentry=tablback[k+1]
                            if firstentry<secondentry
                                tablback[k]=secondentry
                                tablback[k+1]=firstentry
                                sign*=-1
                            else
                                global counter +=1
                            end
                        end
                    end
                    loc=Int64(locfind[tablback])
                    global Q_temp[loc,i]+=copy(sign)
                end
            end
        end
        """transpose in this case is the same as the h.c."""
        Q_temp2=Q_temp#+transpose(Q_temp)
        push!(Qlist,copy(Q_temp2))
    end
end
"""
pervector is a list that I fill with the S matrices as they are created.
But to make it easier to multiply all of them it starts as a list of identities
"""

if Ton=="Hop"
    numQ=length(Qlist)
    Qlist[1]=Ham
    pervector=fill(diagm(ones(ComplexF64,bin)),Int64(numQ))
    pervalue=fill(ones(ComplexF64,bin),Int64(numQ))
    """ diagonalize Qs one by one """
    global S=diagm(ones(ComplexF64,bin))
    for n in 1:Int64(numQ)
        """ diagonalize the product of the previous S matrices with a new Q """
        global decom = eigen(round.(adjoint(S)*Qlist[n]*S, digits=10))#*Q_bunch)#*inv(shiftmat))
        global pervalue[n] = decom.values
        """ save the latest S matrix """
        global pervector[n] = decom.vectors
        """ multiply all the previous s-matrices """
        global S=diagm(ones(ComplexF64,bin))
        for k in 1:Int64(numQ)
            global S*=pervector[k]
        end
    end
    vector=copy(S)
    """
    first we build a set of arrays to save the eigenvalues. every array contains the eigenvalues on a given position.
    """
    combieiglist=fill(ones(Int64(numQ)),bin)
    for i in 1:bin
        global combieig=zeros(Int64(numQ))
        for j in 1:Int64(numQ)
            global combieig[j]=((real.(diag(inv(S)*Qlist[j]*S))))[i]
        end
        combieiglist[i]=combieig
    end

    """
    Oepsies is made to check for degeneracy.
    If 2 sets of eigenvalues are the same the indices of these two will be added to oepsies
    """
    oepsies=[]
    for k in 1:bin-1
        for l in k+1:bin
            same=0
            for m in 1:Int64(numQ)
                if combieiglist[k][m] == combieiglist[l][m]
                    same+=1
                end
            end
            if same == Int64(numQ)
                push!(oepsies,[k,l])
            end
        end
    end
    if length(oepsies)!=0
        print(oepsies)
        error("degeneracies found")
    end

    """
    momdict contains unique sets of M different values between 1 and N
    """
    momdict = collect(combinations(1:N,M))
    for i in 1:bin
        momdict[i]=reverse(momdict[i])
    end

    """
    here the different combinations of momenta are tested. for every rank of Q/J(where J are the current like operators which in the code are the Q_l for l>N/2) we create sum_n cos(rank_q*k_n)/sum_n sin(rank_j*k_n).
    And check of these are the same as the eigenvalues found through diagonalization earlier.
    for rank_q=1 we check the eigenvalues found through Q_1 etc.
    """
    failures=[]
    results=[]
    founds=[]
    global found=1
    for j in 1:bin

        ranky=1
        global found=0
        global list=collect(Int64, 1:1:bin)
        lengthlistmax=bin
        lengthlistmin=0
        # print("ym num",j)
        #println("new cos ", j)
        for ranky in 1:Int64(N/2)
            """ probe is how we temporarily save the sum of the cosines (or sines) """
            global probe=0

            for i in 1:M
                global probe+=cos(ranky*pichart[Int64(momdict[j][i])])
            end
            """tlsit takes only the values that where added to the list last round. This is a bit of a roundabout way of doing it... """
            global tlist=copy(list[lengthlistmin+1:lengthlistmax])
            #println(probe, tlist)
            for k in tlist
                """
                here we check if the sum of the cos is the same as an eigenvalue.
                We save the positions of the eigenvalues where this is true to the list and then check again for the next rank.
                filtering them out one by one
                """

                if round(probe,digits=3)==round(combieiglist[k][ranky],digits=3)/2
                    push!(list,k)
                    # print("eureka")

                end
            end
            lengthlistmin=copy(lengthlistmax)
            lengthlistmax=length(list)
            # print("cos",ranky)
        end
        """same as above but for the sines"""
        #println("new sin ", j)
        for rankys in 1:Int64(N/2-1)
            # print("sin",ranky)
            global probe=0

            for i in 1:M
                global probe+=sin(rankys*pichart[Int64(momdict[j][i])])
            end
            global tlist=copy(list[lengthlistmin+1:lengthlistmax])
            #println(probe, tlist)
            for k in tlist
                if round(probe,digits=3)==round(combieiglist[k][Int64(rankys+N/2)],digits=3)/2
                    push!(list,k)
                    # print("eureka")
                end
            end
            lengthlistmin=copy(lengthlistmax)
            lengthlistmax=length(list)
        end
        #println(copy(list[lengthlistmin+1:lengthlistmax]))
        if length(copy(list[lengthlistmin+1:lengthlistmax]))==0
            push!(failures, j)
        else
            push!(results,[j,copy(list[end])])
            push!(founds,copy(list[end]))
        end

    end
    Fplaces=setdiff(collect(Int64, 1:1:bin),founds)
    """printing failures gives the j that didn't find a partner Links
    printing Fplaces does the opposite
    """



    """add the sum of the momentumindices to the result for ordering"""
    momlist=[]
    for i in 1:length(results)
        tempy=0
        tempo=0
        for j in 1:M
            tempy+=pichart[Int64(momdict[results[i][1]][j])]
            tempo+=momdict[results[i][1]][j]
        end
        push!(results[i],copy(tempo))
        push!(momlist,copy(tempo))
    end



    """
    With this we now know which collumns to exchange to get the correct S.
    This next loop makes a list of exchanges that must be made
    """
    exchanges=[]
    momom=copy(momlist)
    res=copy(results)
    momaxis=[]
    global indexje=1
    while length(exchanges)!=(bin-length(failures))
        #print("floebel",length(momom))
        counti=1
        while counti-1<length(momom)
            a=findall(momom.<momom[counti])
            b=findall(momom.==momom[counti])
            if length(momom)==1
                push!(exchanges,[res[1][2],indexje])
                push!(momaxis,copy(res[1][3]))
                deleteat!(res,1)
                deleteat!(momom,1)
            elseif length(a)==0  && length(b)!=1
                for l in 1:length(b)
                    push!(exchanges,[res[b[l]][2],indexje])
                    push!(momaxis,copy(res[b[l]][3]))
                    global indexje+=1
                end
                deleteat!(res,b)
                deleteat!(momom,b)
            elseif length(a)==0 && length(b)==1
                push!(exchanges,[res[counti][2],indexje])
                push!(momaxis,copy(res[counti][3]))
                global indexje+=1
                deleteat!(res,counti)
                deleteat!(momom,counti)

            end
            counti+=1
        end
    end

    for koppel in exchanges
        vector[:,koppel[2]]=copy(S[:,koppel[1]])
    end

else
    global decom = eigen(Ham)
    global eigvalues=decom.values
    global S=decom.vectors
    vector=copy(S)

end
#
#     combieiglist=fill(ones(Int64(numQ)),bin)
#     for i in 1:bin
#         global combieig=zeros(Int64(numQ))
#         for j in 1:Int64(numQ)
#             global combieig[j]=((real.(diag(inv(S)*Ham*S))))[i]
#         end
#         combieiglist[i]=combieig
#     end
#
#     """
#     Oepsies is made to check for degeneracy.
#     If 2 sets of eigenvalues are the same the indices of these two will be added to oepsies
#     """
#     oepsies=[]
#     for k in 1:bin-1
#         for l in k+1:bin
#             same=0
#             for m in 1:Int64(numQ)
#                 if round(combieiglist[k][m],digits=8) == round(combieiglist[l][m],digits=8)
#                     same+=1
#                 end
#             end
#             if same == Int64(numQ)
#                 push!(oepsies,[k,l])
#             end
#         end
#     end
#     if length(oepsies)!=0
#         print(oepsies)
#         error("degeneracies found")
#     end
# end

vectororder = zeros(ComplexF64,(bin,bin))


if Ton=="hop" && repuls==0 && M==1 && doeik==1
    eiks=Dict()
    eiksum=[]

    for r in ProgressBar(1:bin)
        etemp=0
        for p in 1:M
            # for meh 1:M
            #     change sign
            #println(taart[Int64(ttotal[r][p])])
            etemp+=eiks[taart[Int64(ttotal[r][p])]]
            #println(ttotal[r][p])
        end
        push!(eiksum,etemp/M)
    end

    for r in 1:N
        eiks[taart[r]]=exp(taart[r]im*pi)
    end
    for t in 1:N
        global places=findall(round.(eiksum.*100).==round(moms[t]*100))
        #println(round(moms[t]*100),places)
        global vectororder[:,places[1]] = vector[:,t]
    end
else
    global vectororder=vector
end

println("start prolist")
projs=[string(0,"-",M)]
for p in 1:M
    push!(projs,string(p,"-",M-p))
end

if half==1
    print("huh")
    if aim==0
        global i=1
        happy=0
        protot=zeros(bin,bin)
        global prolist=Dict()
        global spot=round((sizepar+0.0001)/2)
        if allop==1
            print("--------------------------------------------------- hiero")
            for l in 1:length(sitenump)
                global pro=zeros(bin,bin)
                global happen=0
                for k in 1:bin

                    tablcheck=copy(ttotal[k])
                    #if any(spot-sizepar/2 .<= tablcheck .<= spot+sizepar/2))
                    if any(length(findall(spot-sizepar/2 .<= tablcheck .< spot+sizepar/2)) .== sitenump[l]) && (spot-sizepar/2)>0 && (spot+sizepar/2)<=N
                        global pro[k,k]=1
                        global happen=1
                        #print("yes!",tablcheck)
                    elseif (spot-sizepar/2)<0 && (spot+sizepar/2)>=N
                        error("sizepar larger than system")
                    elseif (spot-sizepar/2)<=0 &&  any((length(findall(tablcheck .< spot+sizepar/2))+length(findall(tablcheck .>= N+spot-sizepar/2))) .== sitenump[l])
                        global pro[k,k]=1
                        global happen=1
                    elseif (spot+sizepar/2)>N && any((length(findall(tablcheck .>= spot-sizepar/2))+length(findall(tablcheck .< spot+sizepar/2-N))) .== sitenump[l])
                        global pro[k,k]=1
                        global happen=1
                    end
                end

                if happen!=0
                    global protot+=copy(pro)
                    prolist[l]=copy(pro)
                    # print("\n")
                    # print(i)
                    print("------------------------------------------------------------------------------------------------------happy1")
                end
            end
            if tr(I-protot) != 0
                print("Whoops\n")
                prolist[length(sitenump)+1]=I-protot
            end
        else
            print(" of toch hiero?" )
            global pro=zeros(bin,bin)

            for k in 1:bin

                tablcheck=copy(ttotal[k])
                #if any(spot-sizepar/2 .<= tablcheck .<= spot+sizepar/2))
                if any(length(findall(spot-sizepar/2 .<= tablcheck .< spot+sizepar/2)) .== sitenump) && (spot-sizepar/2)>0 && (spot+sizepar/2)<=N
                    global pro[k,k]=1
                    global happen=1
                    #print("yes!",tablcheck)
                elseif (spot-sizepar/2)<0 && (spot+sizepar/2)>=N
                    error("sizepar larger than system")
                elseif (spot-sizepar/2)<=0 &&  any((length(findall(tablcheck .< spot+sizepar/2))+length(findall(tablcheck .>= N+spot-sizepar/2))) .==sitenump)
                    global pro[k,k]=1
                    global happen=1
                elseif (spot+sizepar/2)>N && any((length(findall(tablcheck .>= spot-sizepar/2))+length(findall(tablcheck .< spot+sizepar/2-N))) .==sitenump)
                    global pro[k,k]=1
                    global happen=1
                end
            end
            if happen!=0
                global protot+=pro
                global prolist[1]=pro
                # print("\n")
                # print(i)
                print("------------------------------------------------------------------------------------------------------happy2")
            end
            if tr(I-protot) != 0
                print("Whoops\n")
                prolist[2]=I-protot
            end
        end

    else
        if average==1
            global bestp=potweight+repuls
            for i in 1:N
                if i+sizepar-1<=N && (sum(pot[i:i+sizepar-1])/sizepar)<bestp
                    global bestp=(sum(pot[i:i+sizepar-1])/sizepar)
                    global spot=floor(sizepar/2)+i
                    # print(i," ",  spot, " ", bestp, "\n")
                elseif i+sizepar-1>N && (sum(pot[i:end])+sum(pot[1:(i+sizepar-N-1)]))/sizepar<bestp
                    global bestp=(sum(pot[i:end])+sum(pot[1:i+sizepar-N-1]))/sizepar
                    global spot=floor(sizepar/2)+i
                    if spot>N
                        spot-=N
                    end
                    # print(i," ",  spot, " ", bestp, "\n")
                end
            end
        else
            spot=argmin(pot)
        end
        global i=1
        happy=0
        protot=zeros(bin,bin)
        prolist=Dict()
        global pro=zeros(bin,bin)

        for k in 1:bin

            tablcheck=copy(ttotal[k])
            #if any(spot-sizepar/2 .<= tablcheck .<= spot+sizepar/2))
            if any(length(findall(spot-sizepar/2 .<= tablcheck .< spot+sizepar/2)) .== sitenump) && (spot-sizepar/2)>0 && (spot+sizepar/2)<=N
                global pro[k,k]=1
                global happen=1
                #print("yes!",tablcheck)
            elseif (spot-sizepar/2)<0 && (spot+sizepar/2)>=N
                error("sizepar larger than system")
            elseif (spot-sizepar/2)<=0 &&  any((length(findall(tablcheck .< spot+sizepar/2))+length(findall(tablcheck .>= N+spot-sizepar/2))) .== sitenump)
                global pro[k,k]=1
                global happen=1
            elseif (spot+sizepar/2)>N && any((length(findall(tablcheck .>= spot-sizepar/2))+length(findall(tablcheck .< spot+sizepar/2-N))) .== sitenump)
                global pro[k,k]=1
                global happen=1
            end
        end
        if happen!=0
            global protot+=pro
            prolist[1]=pro
            print("\n")
            print(i)
            print("------------------------------------------------------------------------------------------------------happy3")
        end

        if tr(I-protot) != 0
            print("Whoops\n")
            prolist[2]=I-protot
        end
    end
else

    if aim==0
        print("hier?")
        global i=1
        happy=0
        protot=zeros(bin,bin)
        global prolist=Dict()
        global spot=round(sizepar/2)
        if allop==1

            for l in 1:length(sitenump)
                global pro=zeros(bin,bin)
                global happen=0
                for k in 1:bin

                    tablcheck=copy(ttotal[k])
                    #if any(spot-sizepar/2 .<= tablcheck .<= spot+sizepar/2))
                    if any(length(findall(tablcheck .% 2 .== 0)) .== sitenump[l])
                        global pro[k,k]=1
                        global happen=1
                        #print("yes!",tablcheck)

                    end
                end

                if happen!=0

                    global protot+=copy(pro)
                    prolist[l]=copy(pro)
                    print("\n")
                    print(i)
                    print("------------------------------------------------------------------------------------------------------happy4")
                end
            end
            if tr(I-protot) != 0
                print("Whoops\n")
                prolist[length(sitenump)+1]=I-protot
            end
        else
            print(" --- of toch hiero? ------" )
            global pro=zeros(bin,bin)

            for k in 1:bin

                tablcheck=copy(ttotal[k])
                #if any(spot-sizepar/2 .<= tablcheck .<= spot+sizepar/2))
                if any(length(findall(tablcheck .% 2 .== 0)) .== sitenump[2:M+1])# && (spot-sizepar/2)>0 && (spot+sizepar/2)<=N
                    global pro[k,k]=1
                    global happen=1
                    #print("yes!",tablcheck)

                end
            end
            if happen!=0
                global protot+=pro
                global prolist[1]=pro
                print("\n")
                print(i)
                print("------------------------------------------------------------------------------------------------------happy5")
            end
            if tr(I-protot) != 0
                print("Whoops\n")
                prolist[2]=I-protot
            end
        end




    else
        print("hier dan?")
        if average==1
            global bestp=potweight+repuls
            for i in 1:N
                if i+sizepar-1<=N && (sum(pot[i:i+sizepar-1])/sizepar)<bestp
                    global bestp=(sum(pot[i:i+sizepar-1])/sizepar)
                    global spot=floor(sizepar/2)+i
                    # print(i," ",  spot, " ", bestp, "\n")
                elseif i+sizepar-1>N && (sum(pot[i:end])+sum(pot[1:(i+sizepar-N-1)]))/sizepar<bestp
                    global bestp=(sum(pot[i:end])+sum(pot[1:i+sizepar-N-1]))/sizepar
                    global spot=floor(sizepar/2)+i
                    if spot>N
                        spot-=N
                    end
                    # print(i," ",  spot, " ", bestp, "\n")
                end
            end
        else
            spot=argmin(pot)
        end
        global i=1
        happy=0
        protot=zeros(bin,bin)
        prolist=Dict()
        global pro=zeros(bin,bin)

        for k in 1:bin

            tablcheck=copy(ttotal[k])
            #if any(spot-sizepar/2 .<= tablcheck .<= spot+sizepar/2))
            if any(length(findall(spot-sizepar/2 .<= tablcheck .< spot+sizepar/2)) .== sitenump) && (spot-sizepar/2)>0 && (spot+sizepar/2)<=N
                global pro[k,k]=1
                global happen=1
                #print("yes!",tablcheck)
            elseif (spot-sizepar/2)<0 && (spot+sizepar/2)>=N
                error("sizepar larger than system")
            elseif (spot-sizepar/2)<=0 &&  any((length(findall(tablcheck .< spot+sizepar/2))+length(findall(tablcheck .>= N+spot-sizepar/2))) .== sitenump)
                global pro[k,k]=1
                global happen=1
            elseif (spot+sizepar/2)>N && any((length(findall(tablcheck .>= spot-sizepar/2))+length(findall(tablcheck .< spot+sizepar/2-N))) .== sitenump)
                global pro[k,k]=1
                global happen=1
            end
        end
        if happen!=0
            global protot+=pro
            prolist[1]=pro
            print("\n")
            print(i)
            print("------------------------------------------------------------------------------------------------------happy6")
        end

        if tr(I-protot) != 0
            print("Whoops\n")
            prolist[2]=I-protot
        end
    end
    println("hier? of zelfs")
end
println("prolist done")
if half == 1
    projectortype = " half"
else
    projectortype = " stagger"
end

""" here the partial overlaps are calculate """
diads=zeros((length(prolist)+1)*(maxorder-minorder+1),bin-offset)
heater=zeros((ending-start,ending-start))
maincont=zeros((ending-start,ending-start))
maxperlevel=zeros((M+1,bin))
N_l=zeros(ComplexF64,(ending-start,ending-start))
hlist=zeros(bin*(bin-1))
Nlist=zeros(bin*(bin-1))
for order in minorder:maxorder

    indexi=1
    for k in ProgressBar(1:ending)
        vector1=vectororder[:,k]
        indexj=1
        #print("\n",k)
        for l in 1:ending
            vector2=vectororder[:,l]


            """
            here all the options for the different ways of calculating d are initialized.
            These are mostly different normalizations. Currently I am using d which is not normalized
            """
            d=0
            global nlterm=0
            d2=0
            # dmok=0
            # dcompa=0
            # dcompb=0
            # dcompcs=0
            # dcomagein=0
            maxcont=0
            """
            this loop calculates the different partial overlaps and adds them up
            """
            #println((k,l))
            for pronum in 1:length(prolist)
                #comp_right=prolist[pronum]*vector2
                paro_right=prolist[pronum]*vector1
                paro_eind=adjoint(vector2)*paro_right
                # compa_eind=(np.transpose(comp_right)).dot(paro_right)
                # compb_eind=(vector2).dot(comp_right)
                #compalta_eind=transpose((paro_right))*paro_right
                #compaltb_eind=transpose((comp_right))*comp_right
                ## print(paro_eind[0,0])
                d += abs(paro_eind)^order
                global nlterm += (pronum-1)*(paro_eind)
                if k==l
                    diads[(order-minorder)*(length(prolist)+1)+pronum,l]= abs(paro_eind)^order
                end
                if abs(paro_eind) > d2
                    d2=abs(paro_eind)
                    maxcont=pronum
                end
                #d2 += ((abs(paro_eind/(compalta_eind*compaltb_eind)^0.5))^order)/length(prolist)
            end
            if k>offset
                diads[(order-minorder)*(length(prolist)+1)+length(prolist)+1,k-offset]=sum(diads[(order-minorder)*(length(prolist)+1)+1:(order-minorder)*(length(prolist)+1)+length(prolist),k-offset])
            end
            """
            when the sum of the partial overlap is completed it gets inserted in a martix which we can later use to plot a graph
            """
            heater[indexi,indexj]=d

            N_l[indexi,indexj]+=nlterm
            if indexi!=indexj
                hlist[Int64((indexi-1)*(bin-1)+(indexj))]=abs(d)
                Nlist[Int64((indexi-1)*(bin-1)+(indexj))]=abs(nlterm)
            end
            maincont[indexi,indexj]=maxcont
            indexj+=1
        end
        indexi+=1
    end

    #if order==1
    global heat=copy(heater)
    #end
    fig = Figure()
    ax,hm = heatmap(fig[1,1], eigvalues, eigvalues, heater, colormap=:hot, title=string(round.(pot*100)/100), aspect=1)
    colsize!(fig[1,1].layout, 1, Aspect(1, 21))
    Colorbar(fig[1,2],hm)
    ax.aspect=1
    ax.yreversed = true
    ax.xlabel="Energy"
    ax.ylabel="Energy"
    supertitle = Label(fig[0, :], string(round.(pot*100)/100), textsize = 20)
    display(fig)
    #save(string("C:\\Users\\reime\\OneDrive\\Desktop\\MasterThesis\\Presentation\\difbw\\energy_scale_",Ton,"_",order,"_rep=",repuls,"_potweight=",potweight,"_", N, "_", M, "_", projectortype,"_pots=", halfpot, "_bw=", bonusweight,"_seed=",currentseed,".pdf"),fig)

    #title!(string(round.(pot*100)/100),fontsize=1)
    #savefig(string("C:\\Users\\reime\\OneDrive\\Desktop\\MasterThesis\\Presentation\\difbw\\",Ton,"_",order,"_rep=",repuls,"_potweight=",potweight,"_", N, "_", M, "_", projectortype,"_pots=", halfpot, "_bw=", bonusweight,"_seed=",currentseed,".pdf"))
    #display(plot(heatmap(heater, yflip=true, xlabel=eigvalues color=:hot, aspect_ratio=1)))
    #display(heatmap(eigvalues, eigvalues, heater, yflip=true, color=:hot, aspect_ratio=1))
    #savefig(string("C:\\Users\\reime\\OneDrive\\Desktop\\MasterThesis\\Half+pot=1_",order,".pdf"))
    if order==2
        f = Figure()
        ax2 = Axis(f[1, 1], xlabel="energy level number", ylabel="partial norm")
        scatter!(ax2, diag(heater))
        f
        #savefig(string("C:\\Users\\reime\\OneDrive\\Desktop\\MasterThesis\\Presentation\\Norm_",order,"_",Ton,"_order_rep=",repuls,"_potweight=",potweight,"_", N, "_", M, "_", projectortype,"_pots=", halfpot, "_bw=", bonusweight,"_seed=",currentseed,".pdf"))
    end
end
for t in 1:bin
    tempheat=copy(heater[:,t])
    tempheat[t]-=1

end
#end


print("plotting")
heaty=copy(heat)

df=DataFrame(vec= Int[], partner = Int[], locs= Array[], amps=Array[] )
numchecks=6
doms=Any[]
cutoff=0.2
for i in 1:bin
    tempdom=Any[]
    global checks=zeros(numchecks)
    for j in 1:numchecks
        tempmax=findmax((heat-I)[:,i])
        global checks[j]=copy(tempmax[2])
        heat[tempmax[2],i]-=tempmax[1]
    end
    for k in checks
        cut=1
        amps=copy(vectororder[:,Int64(k)].*vectororder[:,Int64(i)])
        tempamp=findmax(abs.(amps))[1]
        templocs=[]
        tempamps=[]
        while cut>cutoff
            maxamp=findmax(abs.(amps))
            cut=maxamp[1]/tempamp
            if cut>cutoff
                push!(tempamps,(amps[maxamp[2]]))
                push!(templocs,Int64(maxamp[2]))
            end
            tempamp=maxamp[1]
            amps[maxamp[2]]-=amps[maxamp[2]]
        end
        push!(df, Dict(:vec =>i, :partner=>k, :locs=>templocs, :amps=>tempamps))
    end
end
diahi=findall(diag(heater) .> 0.9)
diamid=findall(0.5 .< diag(heater) .< 0.9)
#diamidlo=findall(0.78 .< diag(heater) .< 0.9)
dialo=findall(diag(heater).<0.5)
"""filter([:vec, :partner] => (vec, partner) -> vec == 252 && partner == 251, df)"""


""" the playground starts here """
# #savefig(string("C:\\Users\\reime\\OneDrive\\Desktop\\MasterThesis\\parts_in_first_half=",sitenump,"_",M,"_particles_with_" ,Ton,"_numsites=",N,"_NN=",repuls,".pdf"))
# a=transpose(abs.(pervector[:,1+1]))*(abs.(reverse(pervector[:,end-1])))
# print(a)
# print("bla")
# for l in 1:M+1
#     listofstates=findall(diag(prolist[l]) .== 1.0)
#     usedstates=[]
#     z=0
#     for k in listofstates
#         z+=1
#         push!(usedstates,reverse(ttotal[k]))
#         print(locfind[reverse(usedstates[i])])
#         print("\n")
#     end
#     print("next")
# end
# # a=[0,0,0,0,0,0,0,0]
# # for i in 1:bin
# #     global a=(a .+ (diads[:,i] ./ diads[:,574]))
# # end
# #
# # print(a)
#
# for i in 1:bin-2
#     if round(pervalue[i]*1000)!=round(pervalue[i+1]*1000) && round(pervalue[i+1]*1000)!=round(pervalue[i+2]*1000)
#        print(i)
#        print("\n")
#     end
# end
#
# averagemax=ones(bin)
# for j in 1:bin
#     global vecy=vectororder[:,j]
#     global dist=zeros(N)
#     for i in 1:bin
#         global dist[ttotal[i]].+=(vecy[i])^2
#     end
#     global tempdist=copy(dist)
#     for k in 1:M
#         global tempymax=findmax(tempdist)
#         global averagemax[j]*=tempymax[1]
#         global tempdist[tempymax[2]]=0
#     end
# end
# plot(bar(dist),xticks=1:1:N,xlabel="site",ylabel="weight")
#
#
# avhi=findall(averagemax .> 0.96)
# avmidhi=findall(0.8 .< averagemax .< 0.96)
# avmidlo=findall(0.68 .< averagemax .< 0.8)
# avlo=findall(averagemax.<0.68)

# averagemax=zeros(bin)
# for j in 1:bin
#     global vecy=vectororder[:,j]
#     global dist=zeros(N)
#     for i in 1:bin
#         global dist[ttotal[i]].+=(abs(vecy[i]))^2
#     end
#     global tempdist=copy(dist)
#     for k in 1:M
#         global tempmax=findmax(tempdist)
#         global averagemax[j]+=tempmax[1]
#         global tempdist[tempmax[2]]=0
#     end
#     averagemax[j]/=M
# end
# plot(bar(dist),xticks=1:1:N,xlabel="site",ylabel="weight")

#
# avhi=findall(averagemax .> 0.96)
# avmidhi=findall(0.8 .< averagemax .< 0.96)
# avmidlo=findall(0.68 .< averagemax .< 0.8)
# avlo=findall(averagemax.<0.68)

cut=copy(diag(heater))
cut[30:end].=0.0

#hotstate=findall(diagm(cut) .> 0.01)
#hotstate=findall(0.512 .< (heaty-I))# .< 0.514)
hotstate=findall(maximum(heaty-I)+0.002 .> (heaty-I) .> maximum(heaty-I)-(2*maximum(heaty-I)/10))
# highs=findall((diag(heater)./averagemax).>1.002)
# midsh=findall(0.99.<(diag(heater)./averagemax).<1.002)
#hotstate=[[98,98]]
pastentry=[[0,0]]
entry=0
disttot=zeros(N)
while entry < length(hotstate)
    println("in while")
    global entry+=1
    while [hotstate[Int64(entry)][2],hotstate[Int64(entry)][1]] in pastentry
        global entry+=1
        if entry > length(hotstate)
            break
        end
    end

    if entry > length(hotstate)
        break
    end

    push!(pastentry, [hotstate[Int64(entry)][1],hotstate[Int64(entry)][2]])



    println(hotstate[Int64(entry)][1]," ", hotstate[Int64(entry)][2])
    vec1=vectororder[:,hotstate[Int64(entry)][1]]
    global dist1=zeros(N)

    for i in 1:bin
        global dist1[ttotal[i]].+=adjoint(vec1[i])*vec1[i]
    end
    disttot.+=dist1

    vec2=vectororder[:,hotstate[Int64(entry)][2]]
    global dist2=zeros(N)

    for j in 1:bin
        global dist2[ttotal[j]].+=adjoint(vec2[j])*vec2[j]
    end
    disttot.+=dist2

    global veco=abs.(vectororder[:,hotstate[Int64(entry)][1]].*vectororder[:,hotstate[Int64(entry)][2]])

    tresh=0.1
    tr=1
    numdis=1

    global countmax=findmax(veco)
    global distweight=[countmax[1]]
    global dist=[countmax[2]]
    veco[countmax[2]]=0.0
    while tr>tresh
        numdis+=1
        println("in second while ",numdis)
        countmax=findmax(veco)
        push!(distweight,countmax[1])
        push!(dist,countmax[2])
        veco[countmax[2]]=0.0
        tr=distweight[numdis]/distweight[1]
    end
    setdiff!(distweight,[distweight[end]])
    setdiff!(dist,[dist[end]])

    global disttest=dist1.*dist2
    global domdists=vec1.*vec2
    parties=zeros(length(prolist))
    for pronum in 1:length(prolist)

        paro_right=prolist[pronum]*vec1
        paro_eind=adjoint(vec2)*paro_right

        parties[pronum]= abs(paro_eind)

    end
    # disttot.+=dist
    #display(plot(bar(dist1,title=string(hotstate[Int64(entry)][1]),xlabel="site",xticks=1:1:N),bar(dist2,title=string(hotstate[Int64(entry)][2]),xlabel="site",xticks=1:1:N),bar(distweight, title=string(hotstate[Int64(entry)][1]," ",hotstate[Int64(entry)][2]),xlabel="distribution", xrotation=90, xticks=(1:length(dist), string.(dist))), bar(parties,title=string("partial overlap (=",round(sum(parties),digits=5),")"),xlabel="projector number",xticks=(1:1:6,projs)),ylabel="weight",legend =false ))
    #savefig(string("C:\\Users\\reime\\OneDrive\\Desktop\\MasterThesis\\Presentation\\overlap states\\",N,"_",M,"\\",halfpot,"\\new\\",half,"\\",halfpot,"_Seed=",currentseed,"_partover=",round(sum(parties),digits=5),"_half=",half,"_",hotstate[Int64(entry)][1],"_",hotstate[Int64(entry)][2],".pdf"))
end
#display(plot(bar(disttot),xticks=1:1:N,xlabel="site", title="tot",ylabel="weight",legend =false ))


""" Operators we might like"""

nj_list=[]
for i in 1:N
    global n_j=zeros(bin,bin)
    for j in 1:bin
        if i in ttotal[j]
            global n_j[j,j]=1

        end
    end
    push!(nj_list,copy(n_j))
end


nj_heat=[]
for i in 1:N
    tempnjheat=zeros(ComplexF64,bin,bin)
    for k in 1:bin
        for j in 1:bin
            tempnjheat[k,j]=(adjoint(vectororder[:,k])*nj_list[i]*vectororder[:,j])
        end
    end
    push!(nj_heat,copy(tempnjheat))
end

Neighcorr_list=[]
for i in 1:N-1
    push!(Neighcorr_list,copy(nj_list[i]*nj_list[i+1]))
end
push!(Neighcorr_list,copy(nj_list[N]*nj_list[1]))

neighc_heat=[]
for i in 1:N
    tempneighcheat=zeros(ComplexF64,bin,bin)
    for k in 1:bin
        for j in 1:bin
            tempneighcheat[k,j]=(adjoint(vectororder[:,k])*Neighcorr_list[i]*vectororder[:,j])
        end
    end
    push!(neighc_heat,copy(tempneighcheat))
end
tempneighctheat=zeros(ComplexF64,bin,bin)
for k in 1:bin
    for j in 1:bin
        tempneighctheat[k,j]=(adjoint(vectororder[:,k])*sum(Neighcorr_list)*vectororder[:,j])
    end
end
push!(neighc_heat,copy(tempneighctheat))


njt_list=[]
njt_heat=[]
tev=(exp((-1im).*Ham))
for i in 1:N
    tempnjtheat=zeros(ComplexF64,bin,bin)
    push!(njt_list,adjoint(tev)*nj_list[i]*tev)
    for k in ProgressBar(1:bin)
        for j in 1:bin
            tempnjtheat[k,j]=(adjoint(vectororder[:,k])*njt_list[i]*nj_list[i]*vectororder[:,j])
        end
    end
    push!(njt_heat,copy(tempnjtheat))
end

heat_graph=reshape(heaty-diagm(diag(heaty)),(bin*bin))
for i in 1:N
    nj_graph=abs.(reshape(nj_heat[i]-diagm(diag(nj_heat[i])),(bin*bin)))
    #display(scatter(nj_graph,heat_graph,legend=false,ylabel="partial overlap",xlabel=string("Filling at site ",i),title=string("partition = ",half," potential = ",halfpot)))
    #savefig(string("C:\\Users\\reime\\OneDrive\\Desktop\\MasterThesis\\Presentation\\overlap states\\",N,"_",M,"\\",halfpot,"\\new\\",half,"\\Filling_Seed=",currentseed,"_pot=",halfpot,"_half=",half,".pdf"))
end

for i in 1:N
    nj_graph=abs.(reshape(neighc_heat[i]-diagm(diag(neighc_heat[i])),(bin*bin)))
    if i == N
        next=1
    else
        next=i+1
    end
    #display(scatter(nj_graph,heat_graph,legend=false,ylabel="partial overlap",xlabel=string("neighbour correlation between sites ",i," and ",next),title=string("partition = ",half," potential = ",halfpot)))
    #savefig(string("C:\\Users\\reime\\OneDrive\\Desktop\\MasterThesis\\Presentation\\overlap states\\",N,"_",M,"\\",halfpot,"\\new\\",half,"\\Neighbour_",i,"_Seed=",currentseed,"_pot=",halfpot,"_half=",half,".pdf"))
end

nj_graph=abs.(reshape(neighc_heat[N+1]-diagm(diag(neighc_heat[N+1])),(bin*bin)))
#display(scatter(nj_graph,heat_graph,legend=false,ylabel="partial overlap",xlabel=string("total neighbour correlation"),title=string("partition = ",half," potential = ",halfpot)))
#savefig(string("C:\\Users\\reime\\OneDrive\\Desktop\\MasterThesis\\Presentation\\overlap states\\",N,"_",M,"\\",halfpot,"\\new\\",half,"\\All_neighc_Seed=",currentseed,"_pot=",halfpot,"_half=",half,".pdf"))

njt_list=[]
njt_heat=[]
tev=(exp((-1im).*Ham))
for i in 1:N
    tempnjtheat=zeros(ComplexF64,bin,bin)
    push!(njt_list,adjoint(tev)*nj_list[i]*tev)
    for k in ProgressBar(1:bin)
        for j in 1:bin
            tempnjtheat[k,j]=(adjoint(vectororder[:,k])*njt_list[i]*nj_list[i]*vectororder[:,j])
        end
    end
    push!(njt_heat,copy(tempnjtheat))
end


nj_graph=abs.(reshape(sum(njt_heat[1:Int64(N/2)])-diagm(diag(sum(njt_heat[1:Int64(N/2)]))),(bin*bin)))
#display(scatter(nj_graph,heat_graph,legend=false,ylabel="partial overlap",xlabel=string("time-correlated filling of left half"),title=string("partition = ",half," potential = ",halfpot)))
#savefig(string("C:\\Users\\reime\\OneDrive\\Desktop\\MasterThesis\\Presentation\\overlap states\\",N,"_",M,"\\",halfpot,"\\new\\",half,"\\Left_time_filling_Seed=",currentseed,"_pot=",halfpot,"_half=",half,".pdf"))
nj_graph=abs.(reshape(sum(nj_heat[1:Int64(N/2)])-diagm(diag(sum(nj_heat[1:Int64(N/2)]))),(bin*bin)))
#display(scatter(nj_graph,heat_graph,legend=false,ylabel="partial overlap",xlabel=string("filling of left half"),title=string("partition = ",half," potential = ",halfpot)))
#savefig(string("C:\\Users\\reime\\OneDrive\\Desktop\\MasterThesis\\Presentation\\overlap states\\",N,"_",M,"\\",halfpot,"\\new\\",half,"\\Left_no-time_filling_Seed=",currentseed,"_pot=",halfpot,"_half=",half,".pdf"))

# tempnjtheat=zeros(ComplexF64,bin,bin)
# for k in ProgressBar(1:bin)
#     for j in 1:bin
#         tempnjtheat[k,j]=(adjoint(vectororder[:,k])*sum(njt_list[1:Int64(N/2)])*sum(nj_list[1:Int64(N/2)])*vectororder[:,j])
#     end
# end
#
# nj_graph=abs.(reshape(tempnjtheat-diagm(diag(tempnjtheat)),(bin*bin)))
# display(scatter(nj_graph,heat_graph,legend=false,ylabel="partial overlap",xlabel=string("time-correlated filling of left half"),title=string("partition = ",half," potential = ",halfpot)))
# #savefig(string("C:\\Users\\reime\\OneDrive\\Desktop\\MasterThesis\\Presentation\\overlap states\\",N,"_",M,"\\",halfpot,"\\new\\",half,"\\Left_sum-first-time_filling_Seed=",currentseed,"_pot=",halfpot,"_half=",half,".pdf"))

# jneigh=1
# Neigh=zeros(bin,bin)
# for i in 1:bin
#     for k in 1:bin
#         for l in 1:M
#             p=ttotal[i][l]
#             if p+jneigh > N
#                 p-=N
#             end
#             if p+jneigh in ttotal[k]
#                 Neigh[i,k]+=1
#             end
#         end
#     end
# end


#
# ki=0
# for i in pastentry
#     ki+=1
#     kj=0
#     for j in pastentry
#         kj+=1
#         if i==reverse(j)
#             println(ki," ",kj)
#         end
#     end
# end
dia2=[]
for i in 1:bin
    parties=zeros(M+1)
    for pronum in 1:length(prolist)

        paro_right=prolist[pronum]*vectororder[:,i]
        paro_eind=adjoint(vectororder[:,i])*paro_right

        parties[pronum]= abs(paro_eind)^2

    end
    push!(dia2,sum(parties))
end


vnum=1
vec=vectororder[:,vnum]
dista=zeros(ComplexF64,N)

for i in 1:bin
    dista[ttotal[i]].+=(vec[i])^2
end
#display(plot(bar(abs.(dista)),xticks=1:1:N,xlabel="site", title=string(vnum),ylabel="weight",legend =false ))

vnumb=5
vecb=vectororder[:,vnumb]
distb=zeros(ComplexF64,N)

for i in 1:bin
    distb[ttotal[i]].+=(vecb[i])^2
end
#display(plot(bar(abs.(distb)),xticks=1:1:N,xlabel="site", title=string(vnumb),ylabel="weight",legend =false ))
distel=dista.*distb

grouping=findall(0.6 .< diag(heater) .< 1.1)
distht=zeros(Float64,N)
tik=0
for k in grouping
    global tik+=1
    vech=vectororder[:,k]
    global disth=zeros(Float64,N)
    for i in 1:bin
        global disth[ttotal[i]].+=(vech[i])^2
    end
    #display(bar(disth))
    disth.=real.(disth)
    if tik <= length(grouping)/2
        global distht.+=disth
    elseif tik == length(grouping)/2 + 1
        global tiky=tik
    end
end
distht./=length(grouping)/2
# vec=vectororder[:,1]
# dist1=zeros(N)
# dist2=zeros(N)
#
# for j in 1:bin
#     vec=vectororder[:,j]
#     if j in highs
#         for i in 1:bin
#             dist1[ttotal[i]].+=(vec[i])^2
#         end
#     else
#         for i in 1:bin
#             dist2[ttotal[i]].+=(vec[i])^2
#         end
#     end
# end
# plot(bar(dist1./length(highs)),xticks=1:1:N,xlabel="site",ylabel="weight")
# plot(bar(dist2./length(lows)),xticks=1:1:N,xlabel="site",ylabel="weight")

strength=zeros(bin)
averageoverlap=zeros(bin)
for i in 1:bin
    strength[i]=sum(abs.(vectororder[:,i]))
    averageoverlap[i]=sum(heater[:,i])-1
end



# h=copy(heater)
# weights=zeros(N)
# for i in filter([:vec, :partner] => (vec, partner) -> vec == 182 && partner == 81, df)[!, "locs"][1]
#     weights[ttotal[i]].+=1
# end
