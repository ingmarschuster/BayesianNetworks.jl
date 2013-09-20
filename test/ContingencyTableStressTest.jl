module ContingencyTableStressTest

using DataFrames
using ContingencyTables
using QuickCheck

  function random_table(nvars::Int64)
    init = false
    tab = ContingencyTable(StickBreakingProbabilities(2))
    
    for i = 1:nvars
	data = DataFrame()
	data["probab"]=StickBreakingProbabilities(2)
	data[string("i_",i)] = [0:length(data["probab"])-1]
	smalltab = ContingencyTable(data)
	
	if !init
	    tab = deepcopy(smalltab)
	    init = true
	else
	    tab = fproduct(tab,smalltab)
	end
    end
 
    return tab
  end

  function random_table_create(maxsize::Int64, eachamount::Int64 = 100)
    avtime = Array(Float64,maxsize)
    avtime2 = Array(Float64, maxsize)
    maxdim = Array(Float64,maxsize)
    all_tables = Array(ContingencyTable,0)
    rest = 0
    for j = 1:maxsize
	avtime[j] = 0
	avtime2[j] = 0
	maxdim[j] = 0
	rest = (j-1)*eachamount
	for k = 1:eachamount
	    tic()
	    rand = random_table(j)
	    lines = length(rand.data["probab"])
	    rows = length(colnames(rand.data))
	    println("created ",lines," x ",rows, " ContingencyTable")
	    temp = toc()
	    println("progress for ",j," variables: ",round(k/eachamount*100,2),"%\t(overall ",round((k+rest)/(maxsize*eachamount)*100,2),"%)")
	    push!(all_tables,rand)
	    avtime[j]+=temp
	    avtime2[j]+=temp/(lines*rows)
	    if lines*rows > maxdim[j]
		maxdim[j] = lines*rows
	    end
	    #println("time per field: ", (temp)/(lines*rows))
	end
    end
    
    println("\n********\t\tResults\t\t    *******")
    println("***\tcreated ", eachamount," random ContingencyTables\t***")
    println("***\tfor each 1 to ", maxsize, " random variables\t***") 
    
    for i = 1:maxsize
	println()
	println("average time per ConTable for ",i, " random variables: ", avtime[i]/eachamount*10^3, " ms")
	println("average time per data field for ",i, " random variables: ", avtime2[i]/eachamount*10^3, " ms")
	println("average number of fields per ConTable for ",i," random variables: ",int(round(avtime[i]/avtime2[i]))," (max: ",int(maxdim[i]),")")
	end
    return all_tables
  end
  
  function create_testing_data(maxsize::Int64, eachamount::Int64 = 100)
      file_stream = open("random_data.txt","w")
      serialize(file_stream,random_table_create(maxsize,eachamount))
      close(file_stream)
  end
  
  function stresstest(ntests::Int64)
      file_stream = open("random_data.txt")
      data = deserialize(file_stream)
      close(file_stream)
      
      possible_ev::Array{Function}
      possible_ev= [x->mod(x,2) == 0, x->x==1, x->contains([1,2,3],x)]
      
      for i = 1:ntests
	      x = rand(1:length(data))
	  
	      table1 = data[x]
	      tic()
	      println("processing ",length(table1.data["probab"])," x ",length(colnames(table1.data))," ContingencyTable...")
	      println("conditioning...")
	      condition(table1, [setdiff(colnames(table1.data),["probab"])[rand(1:length(colnames(table1.data))-1)]])
	      println("done")
	      println("marginalizing...")
	      marginalize(table1, setdiff(colnames(table1.data),["probab"])[rand(1:length(colnames(table1.data))-1)])
	      println("done")
	      println("observing evidence...")
	      evidence(table1, [(setdiff(colnames(table1.data),["probab"])[rand(1:length(colnames(table1.data))-1)],possible_ev[rand(1:length(possible_ev))])])
	      println("done")
	      toc()
	      println(ntests-i," tests remaining")
      end
  end
  
  #create_testing_data(5,10000)
  stresstest(100)
end