module ContingencyTables

    using DataFrames
    using Distributions
    import Base.isequal
    
    export JointProbability
    export ContingencyTable
    export fproduct
    export StickBreakingProbabilities
    export condition
    export evidence
    export marginalize
    export normalize!
    const EPSILON = 1/10^10
    
    abstract JointProbability

    type ContingencyTable <: JointProbability
        data::DataFrame
        evidence::Array{(String, Function)}
        conditionedon::Array{String}
        function ContingencyTable(probabilities::Array{Float64},ev = Array((String,Function),0), con = Array(String,0))
            if abs(1-sum(probabilities))>EPSILON && isempty(con)
		#println(abs(sum(probabilities)-1))
		warn("probabilities don't sum to 1")
            end
            df = DataFrame(i = 1:length(probabilities), probab = probabilities)
            new(df,ev,con)
        end
        function ContingencyTable(a::DataFrame,ev = Array((String,Function),0), con = Array(String,0))
	    if abs(1-sum(a["probab"]))>EPSILON && isempty(con)
		#println(abs(sum(a["probab"])-1)
		warn("probabilities don't sum to 1")
	    end
	    new(a,ev,con)
        end
    end

    function StickBreakingProbabilities(concentration::Number = 2)    
        stick = 1.0
        events = rand(Distributions.Poisson(concentration))+2 # there have to be at least two discrete values; else there would be only one RV level with probability 1.
        
    
        retval = Array(Float64, events)
        for i = 1:(events - 1)
            retval[i] = stick * rand() #rand() returns random draws between 1 and 0
            stick -= retval[i]
        end
        retval[events] = stick
        retval
    end
    
    function evidence(table::ContingencyTable, evidence, normalizing)
	#normalizing = true
	if typeof(evidence) != Array{(ASCIIString,Function),1} && typeof(evidence) != Array{(Any...,),1}
	    error("evidence must be of type Array{(ASCIIString, Function)}")
	end
	
	conditioned = DataFrame()
	columns = colnames(table.data)
	
	for pair in evidence
	    columns = setdiff(columns,[pair[1]])
	end
	
	if isequal(columns, colnames(table.data))
	    warn("observed variable doesn't exist")
	    return table
	end
	
	initialized = false
	include = false
	
	for row = 1:nrow(table.data)
	    include = true
	    for pair in evidence
		if contains(colnames(table.data),pair[1])
		    if !pair[2](table.data[row, pair[1]])
			include = false
		    end
		end
	    end

	    if include
		if !initialized
		      conditioned = deepcopy(table.data[row,columns])
		      initialized = true
		else
		      conditioned = rbind(conditioned, table.data[row,columns])
		end
	    end
	end
	
	#conditioned["probab"] = map(x->x/sum(conditioned["probab"]), conditioned["probab"])
	new_tab = ContingencyTable(conditioned, deepcopy(table.evidence), deepcopy(table.conditionedon))
	if normalizing
	  normalize!(new_tab)
	end
	
	for pair in evidence
	    if !contains(new_tab.evidence, pair)
		push!(new_tab.evidence, pair)
	    end
	    
	end
	
	return new_tab
    end
    
    function evidence(table::ContingencyTable, evidence)
	return evidence(table::ContingencyTable, evidence, true)
    end
    
    function evidence!(table::ContingencyTable, evidence)
	table = deepcopy(evidence(table, evidence))
    end
    
    function condition(table::ContingencyTable, conditionedon)
	if typeof(conditionedon) != Array{ASCIIString,1}
	    error("conditionedon must be of type Array{SCIIString}")
	end
	
	if intersect(colnames(table.data), conditionedon) == []
	    warn("condition variable doesn't exist")
	    return table
	end
	
	not_cond = setdiff(setdiff(colnames(table.data), union(conditionedon,["probab"])),table.conditionedon)
	marginalized = deepcopy(table)
	
	for rowname in not_cond
	    marginalized = marginalize(marginalized, rowname)
	end

	condition_vars = setdiff(colnames(marginalized.data), ["probab"])
	conditioned = deepcopy(table)
	
	for i = 1:length(table.data[1])
	    for j = 1:length(marginalized.data[1])
		if isequal(conditioned.data[condition_vars][i,:],marginalized.data[condition_vars][j,:])
		#if sum(map(x->abs(x),values(conditioned.data[condition_vars][i,:])-values(marginalized.data[condition_vars][j,:])))==[0]
		    conditioned.data["probab"][i] = conditioned.data["probab"][i]/marginalized.data["probab"][j]
		end
	    end
	end
	
	for el in conditionedon
	    push!(conditioned.conditionedon, el)
	end
	
	return conditioned
    end
    
    function condition!(table::ContingencyTable, conditionedon)
	table = deepcopy(condition(table,conditionedon))
    end
    
    function marginalize(table::ContingencyTable, key::String)
	if intersect(colnames(table.data), [key]) == []
	    warn("key variable doesn't exist")
	    return table
	end
	
	new_keys = setdiff(colnames(table.data), [key, "probab"])
	include_list = []
	exclude_list = []
	
	for i = 1:length(table.data[1])
	    if !contains(exclude_list, i)
		include_list = vcat(include_list, i)
		for j = i+1:length(table.data[1])
		    if isequal(table.data[new_keys][i,:],table.data[new_keys][j,:])
		    #if sum(map(x->abs(x),values(table.data[new_keys][i,:])-values(table.data[new_keys][j,:])))==[0]
			exclude_list = vcat(exclude_list, j)
		    end
		end
	    end
	end
	
	marginalized = deepcopy(table.data[union(new_keys, ["probab"])][include_list,:])
	rest = 0
	for i = 1:length(include_list)
	    marginalized["probab"][i] = table.data["probab"][include_list[i]]
	    for j = 1:length(exclude_list)/length(include_list)
		marginalized["probab"][i]+= table.data["probab"][exclude_list[j+rest]]
	    end
	    rest += length(exclude_list)/length(include_list)
	end
	
	new_tab = ContingencyTable(marginalized,deepcopy(table.evidence),deepcopy(table.conditionedon))
	
	if intersect(new_tab.conditionedon,[key]) != []
	    new_tab.conditionedon = setdiff(new_tab.conditionedon, [key])
	    for i = 1:length(include_list)
		marginalized["probab"][i] /=2
	    end
	end
	
	return new_tab
    end
    
    function marginalize!(table::ContingencyTable, key::String)
	table = marginalize(table,key)
    end
    
    function fproduct(table1::ContingencyTable, table2::ContingencyTable, normalizing::Bool = true)
 	subtable1 = true
 	for elem in colnames(table2.data)
 	    if !contains(colnames(table1.data),elem)
 		subtable1 = false
 	    end
 	end
 	
 	subtable2 = true
 	for elem in colnames(table1.data)
 	    if !contains(colnames(table2.data),elem)
 		subtable2 = false
 	    end
 	end
	
	product = DataFrame()
	if (subtable1 || subtable2) && isempty(table1.conditionedon) && isempty(table2.conditionedon) && colnames(table1.data) != ["probab"] && colnames(table2.data) != ["probab"]
	    error("ContingencyTables must contain at least 1 different variable")
	else
	    cols = setdiff(intersect(colnames(table1.data),colnames(table2.data)),["probab"])
	    table1_cols = setdiff(colnames(table1.data),["probab"])
	    table2_cols = setdiff(colnames(table2.data), cols)
	    current_column = []
	    pos = 1
	
	    for i = 1:nrow(table1.data)
		for j = 1:nrow(table2.data)
		    #if sum(map(x->abs(x),values(table1.data[i,cols])-values(table2.data[j,cols])))==[0] || isequal(cols,[])
		    if isequal(table1.data[i,cols],table2.data[j,cols]) || isequal(cols,[])
			current_column = hcat(table1.data[i,table1_cols], table2.data[j,table2_cols])
			if isempty(product)
			    product = deepcopy(current_column)
			else
			    product = rbind(product, current_column)
			    pos += 1
			end
			product["probab"][pos] *= table1.data[i, "probab"] 
		    end
		end
	    end
	end
	
	#println(sum(product["probab"]))
	ev = Array((String, Function),0)
	if(table1.evidence!=[] || table2.evidence!=[])
	    for el in table1.evidence
		if !contains(colnames(product),el[1])&&!contains(ev,el)
		    push!(ev,el)
		end
	    end
	    for el in table2.evidence
		if !contains(colnames(product),el[1])&&!contains(ev,el)
		    push!(ev,el)
		end
	    end
	end
	
	product = ContingencyTable(product,ev, unique(union(table1.conditionedon, table2.conditionedon)))
	if normalizing normalize!(product) end
	
	return product
    end
    
    function fproduct!(table1::ContingencyTable, table2::ContingencyTable, normalizing::Bool = true)
	table1 = fproduct(table1,table2, normalizing)
	return true
    end
    
    function normalize!(table::ContingencyTable)
	if isempty(table.conditionedon)
	    table.data["probab"] = map(x->x/sum(table.data["probab"]), table.data["probab"])
	else
	   #partitioning table and normalizing each subtable
	   partition = Array(Array,0)
	   values = Array(DataFrame,0)
	   
	   for i = 1:nrow(table.data)
	      row = table.data[i,table.conditionedon]
	      
	      if contains(values, row)
		  for j = 1:length(values)
		      if isequal(values[j], row)
			  push!(partition[j],i)
		      end
		  end
	      else
		  push!(values, row)
		  push!(partition, [i])
	      end
	   end
	   
	   for array in partition
	      table.data[array,"probab"] = map(x->x/sum(table.data[array,"probab"]), table.data[array,"probab"])
	   end
	end
	return true
    end
    
    function normalize(table::ContingencyTable)
	tab = deepcopy(table)
	normalize!(tab)
	return tab
    end
    
    function isequal(t1::ContingencyTable, t2::ContingencyTable, fuzzy::Bool = false)
	if sort(t1.conditionedon) != sort(t2.conditionedon)
	    return false
	end
	if length(t1.evidence) != length(t2.evidence)
	    return false
	else
	    for i = 1:length(t1.evidence)	   
		if sort(t1.evidence)[i][1] != sort(t2.evidence)[i][1]
		    return false
		end
		if (sort(t1.evidence)[i][2].code.ast.args[[1,2]] != sort(t2.evidence)[i][2].code.ast.args[[1,2]]) || (sort(t1.evidence)[i][2].code.ast.args[3].args[2] != sort(t2.evidence)[i][2].code.ast.args[3].args[2])
		    return false
		end
	    end
	end
	if !isequal(sort(colnames(t1.data)), sort(colnames(t2.data)))    
	    return false
	else
	    for j= 1:length(t1.data["probab"])
		#println(t1.data[j,:])
		if !findmatch(t1.data[j,:],t2.data,fuzzy)
		    return false
		end
	    end
	end
      
	return true
    end
  
    function findmatch(pattern::DataFrame, data::DataFrame,fuzzy::Bool)
	if(isempty(colnames(pattern)))
	    return true
	end
	name = colnames(pattern)[1]
	value = pattern[name][1]
	matched = false
	nrow = 0
      
	for elem in data[name]
	    nrow +=1
	    if name!="probab" || !fuzzy
		if elem == value
		    matched = findmatch(pattern[setdiff(colnames(pattern),[name])], data[setdiff(colnames(pattern),[name])][nrow,:],fuzzy)
		end
	    else
		if abs(elem-value) <= EPSILON
		    matched = findmatch(pattern[setdiff(colnames(pattern),[name])], data[setdiff(colnames(pattern),[name])][nrow,:],fuzzy)
		end
	    end
	    if matched return true
	    end
	end
	return false
    end
end