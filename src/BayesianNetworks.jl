module BayesianNetworks
    using ContingencyTables
    using DataFrames
    using Graphs
    
    import Base.getindex
    
    type BayesianNetwork
	net::AbstractGraph#(String,ContingencyTable),Edge}
	function BayesianNetwork(probs)
	    
	    net = try graph(KeyVertex{(String,JointProbability)},Edge{KeyVertex{(String,JointProbability)}})
	    catch e
		if isa(e,BoundsError)
		    warn("input needs to be an Array")
		    new(graph(KeyVertex{(String,typeof(probs))},Edge{KeyVertex{(String,typeof(probs))}}))
		end
	    end
	    
	    nodes = Array((String, JointProbability, Array{String}),0)
	    
	    for elem in probs
		name = setdiff(colnames(elem.data),union(elem.conditionedon,["probab"]))
		if length(name)!=1 error("probability table can only have one unconditional variable") end
		push!(nodes, (name[1], elem, elem.conditionedon))
	    end	    
	    
	    map((x,y)->add_vertex!(net,KeyVertex{(String,JointProbability)}(y,(x[1],x[2]))),nodes,[1:length(nodes)])
	    
	    vs = vertices(net)
	    edges = 1
	    
	    for i = 1:length(nodes)
		for j = 1:length(nodes)
		    if contains(nodes[i][3],nodes[j][1])
			add_edge!(net,Edge(edges,vs[j],vs[i]))
			edges+=1
		    end
		end
	    end
	    
	    new(net)
	end
	
	function BayesianNetwork(net::AbstractGraph)
	    new(net)
	end
    end
    
    function marginal(bnet::BayesianNetwork, var::String, observed_evidence = Array((ASCIIString, Function),0))
 	if !isempty(observed_evidence)
 	    tables = Array(JointProbability,0)
 	    map(x->push!(tables, evidence(x.key[2], observed_evidence, false)), bnet.net.vertices)
 	    #println(tables)
 	    props = joinedprob(tables)
 	else
 	    props = joinedprob(bnet)
 	end
	#props = joinedprob(bnet)
	
	for name in setdiff(varnames(bnet),[var])
	    props = marginalize(props, name)
	end
	props.conditionedon = Array(String,0)
	normalize!(props)
	return props
    end
    
    function joinedprob(bnet::BayesianNetwork)
	init = deepcopy(bnet.net.vertices[1].key[2])
	
 	if length(bnet.net.vertices) > 1
	    for i = 2:length(bnet.net.vertices)
		init = fproduct(init, bnet.net.vertices[i].key[2],false)
	    end
 	end
	init.conditionedon = Array(String,0)
	return init
    end
    
    function joinedprob(tables::Array{JointProbability})
	init = deepcopy(tables[1])
	if length(tables) > 1
	    for i = 2:length(tables)
		init = fproduct(init, tables[i],false)
	    end
	end
	init.conditionedon = Array(String,0)
	return init
    end
    
    function varelim(bnet::BayesianNetwork, var::String, observed_evidence = Array((ASCIIString, Function),0))
	list = Array((Array{String},JointProbability,Array{String}),0)
	varlist = Array(String,0)
	relvars = prunevars(bnet, var, observed_evidence, false)
	println("relevant variables: ",relvars)
	
	dummy=(["dummy"],ContingencyTable(StickBreakingProbabilities(2)))
	for node in bnet.net.vertices
	    if contains(relvars, node.key[1])
		if node.key[1]!= var
		    next = evidence(node.key[2], observed_evidence,false)
		    
		    for parent in parents(bnet, node)
			if !contains(relvars, parent.key[1])
			    next = marginalize(next,parent.key[1])
			end
		    end
		    
		    push!(list, ([node.key[1]],next,Array(String,0)))
		    push!(varlist, node.key[1])
		else
		    dummy = ([node.key[1]],evidence(node.key[2],observed_evidence,false),Array(String,0))
		end
	    end
	end
	push!(list,dummy)
	push!(varlist, dummy[1][1])
	induced = inducedgraph(bnet)
	
	stop = length(varlist)
	for m = 1:stop
	    println(induced)
	    if m == stop
		tomarginalize = var
	    else 
		#println(pickvar(bnet,list,setdiff(varlist,[var]),induced))
		(tomarginalize,induced) = pickvar(bnet,list,setdiff(varlist,[var]),induced)
		varlist = setdiff(varlist, [tomarginalize])
	    end
	    #println(tomarginalize,"\t",varlist)
	    
	    pair = list[1]
	    index = Array(Int64,0)
	    for i = 1:length(list)
		#println(i,"\t",list[i][1])
		if contains(list[i][1],tomarginalize)
		    pair = list[i]
		    push!(index, i)
		end
	    end
	    
	    vars = [pair[1]]
	    #println(pair[1][1])
	    children = childs(bnet, pair[1][1])
	    npair = pair[2]
	    hiddenvars = pair[3]
	    
	    for child in children
		#println("checking for child ", child.key[1])
		for i = 1:length(list)
		    if !contains(index,i) && !isempty(index)
			#println(list[i][1],"\t",list[i][3])
			if contains(list[i][1],child.key[1]) 
			    npair = fproduct(npair, list[i][2],false)
			    #println("product with ",child.key[1])
			    #println(evidence(list[i][2],observed_evidence,false))
			    push!(index, i)
			    #println(list[i][1],"\t",child.key[1])
			    push!(vars, child.key[1])
			    hiddenvars = vcat(hiddenvars, list[i][3])
			else
			    if contains(list[i][3],child.key[1])
				npair = fproduct(npair, list[i][2],false)
				#println("product with ",child.key[1])
			        push!(index, i)
			        vars = vcat(vars, list[i][1])
			        hiddenvars = vcat(hiddenvars, list[i][3])
			    end
			end
		    end
		end
	    end
	    
	    if tomarginalize == var
		npair.conditionedon = Array(String,0)
		normalize!(npair)
		return npair
	    end
	    
	    #println(length(list))
	    #println(vars)
	    println("maginalizing ",tomarginalize)
	    npair = marginalize(npair, tomarginalize)
	    #println("result: ", npair)
	    
	    #println(index)
	    #println("before: ", list)
	    index = sort(index, rev = true)
	    for j in index
		splice!(list,j)
	    end
	    #println("after: ", list)
	    
	    if length(vars)>1
		vars = setdiff(vars, [tomarginalize])
		push!(list,(vars,npair,vcat(hiddenvars,[tomarginalize])))
		#push!(list,(vars[2:end],npair,vcat(hiddenvars,[vars[1]])))
	    else
		push!(list,(vars,npair,pair[3]))
	    end
	end
	#println(length(list))
	error("target variable ", var," wasn't encountered")
    end
    
    function inducedgraph(bnet::BayesianNetwork)
	induced = graph(KeyVertex{String},Edge{KeyVertex{String}})#,is_directed = false)
	map(x->add_vertex!(induced,KeyVertex{String}(x.index,x.key[1])),bnet.net.vertices)
	map(x->add_edge!(induced, Edge(x.index,induced.vertices[x.source.index],induced.vertices[x.target.index])), bnet.net.edges)
	for node in induced.vertices
	   par = parents(induced, node)
	   #println(node,"\t",par)
	   #println(childs(induced,node))
	   if length(par)>1
	      for par1 in par
		  for par2 in par
		     if par1.key != par2.key && !contains(childs(induced,par1),par2) && !contains(childs(induced,par2),par1)
			  e = Edge(length(induced.edges)+1,par1,par2)			  
			  add_edge!(induced, e)
		     end
		  end
	      end
	   end
	end
	induced.is_directed = false
	return induced
    end
   
    function pickvar(bnet::BayesianNetwork, list, varlist, induced::AbstractGraph)
	for possible_pick in varlist
	    if length(childs(bnet, possible_pick)) == 0
		return possible_pick,induced
	    end
	end
	bestcost = Inf
	bestpick = "dummy"
	besttemp = nothing
	for possible_pick in varlist
	    cost = 0
	    children = childs(bnet, possible_pick)
	    indices = Array(String,0)
	    temp = deepcopy(induced)
	    
	    for i = 1:length(list)
		if contains(list[i][1],possible_pick)
		    indices = union(indices, setdiff(list[i][1],possible_pick))
		end
	    end
	    for child in children
		for j = 1:length(list)
		    if contains(list[j][1],child.key[1])
			indices = union(indices, list[j][1])
		    end
		end
	    end
	    
	    index1 = 0
	    index2 = 0
	    for k = 1:length(indices)
		for l = k+1:length(indices)
		    for m = 1:length(induced.vertices)
			#println(indices[k],"\t",induced.vertices[m].key)
			if indices[k] == induced.vertices[m].key
			    index1 = m
			end
			if indices[l] == induced.vertices[m].key
			    index2 = m
			end
		    end
		    
		    if !contains(childs(induced, induced.vertices[index1]),induced.vertices[index2])
			cost+=1
			add_edge!(temp, Edge(length(temp.edges)+1,induced.vertices[index1],induced.vertices[index2]))
		    end
		end
	    end
	    
	    if cost < bestcost
		bestpick = possible_pick
		bestcost = cost
		besttemp = temp
	    end
	end
	
	return bestpick,besttemp
    end
    
    function prunevars(bnet::BayesianNetwork, target::String, observed_evidence, ischild::Bool, seenvars = Array(String,0))
	activevars = Array(String,0)
	if !ischild
	    seenvars = push!(seenvars, target)
	end
	targetnode = bnet.net.vertices[1]
	
	for node in bnet.net.vertices
	    if node.key[1] == target
		targetnode = node
	    end
	end
	
	has_evidence = false
	for el in observed_evidence
	    if targetnode.key[1] == el[1]
		has_evidence = true
	    end
	end
	
	if !has_evidence
	    if ischild
		activechilds = Array(String,0)
		for child in childs(bnet, targetnode)
		    if !contains(seenvars, child.key[1])
			activechilds = vcat(activechilds, prunevars(bnet,child.key[1],observed_evidence,true,seenvars))
			#println(child.key[1]," grants ", activechilds, " for ", target)
		    end
		end
		#println("total of vars granted for ", target , " are ", activechilds)
		if !isempty(activechilds)
		    push!(activevars, target)
		    activevars = vcat(activevars,activechilds)
		end
	    else
		push!(activevars,target)
		for child in childs(bnet, targetnode)
		    if !contains(seenvars, child.key[1])
			activevars = vcat(activevars, prunevars(bnet,child.key[1],observed_evidence,true,seenvars))
		    end
		end
		for parent in parents(bnet, targetnode)
		    if !contains(seenvars, parent.key[1])
			activevars = vcat(activevars, prunevars(bnet, parent.key[1], observed_evidence, false,seenvars))
		    end
		end
	    end
	else 
	    push!(activevars,target)
	    if ischild
		for parent in parents(bnet, targetnode)
		    if !contains(seenvars, parent.key[1])
			activevars = vcat(activevars, prunevars(bnet, parent.key[1], observed_evidence, false,seenvars))
		    end
		end
	    end
	end
	
	#println("after processing ", target, " as ", ischild, " active vars are ", activevars)
	return unique(activevars)
    end
    
    function varnames(bnet::BayesianNetwork)
	vars = Array(String,0)
	map(x->push!(vars, x.key[1]),bnet.net.vertices)
	return vars
    end
    
    function childs(bnet::BayesianNetwork, varname::String)
	for index = 1:length(bnet.net.vertices)
	    if bnet.net.vertices[index].key[1] == varname
		return out_neighbors(bnet.net.vertices[index],bnet.net)
	    end
	end
	error("variable ", varname, " not found")
    end
    
    function childs(bnet::BayesianNetwork, var::KeyVertex{(String,JointProbability)})
	return childs(bnet, var.key[1])
    end
    
    function childs(g::AbstractGraph, var::KeyVertex{String})
	return out_neighbors(var,g)
    end
    
    function parents(g::AbstractGraph, var::KeyVertex{String})
	list = Array(KeyVertex{String},0)
	for node in g.vertices
	    if contains(out_neighbors(node,g),var)
		push!(list,node)
	    end
	end
	return list
    end
    
    function parents(bnet::BayesianNetwork, var::KeyVertex{(String,JointProbability)})
	list = Array(typeof(bnet.net.vertices[1]), 0)
	for index = 1:length(bnet.net.vertices)
	    if contains(childs(bnet, bnet.net.vertices[index]), var)
		push!(list, bnet.net.vertices[index])
	    end
	end
	return list
    end
    
    function parents(bnet::BayesianNetwork, varname::String)
	pair = KeyVertex{(String, ContingencyTable)}
	for index = 1:length(bnet.net.vertices)
	    if bnet.net.vertices[index].key[1] == varname
		pair = bnet.net.vertices[index]
	    end
	end
	
	return parents(bnet, pair)
    end
    
    function getindex(net::BayesianNetwork, varname::String)
	for i = 1:length(net.content)
 	    if varname == net.content[i].name
 		return net.content[i].value
 	    end
 	end
 	error("variable ", varname, " not found")
    end
     
    function getindex(net::BayesianNetwork, index::Int64)
	return net.content[index].value 
    end
    
    a = BayesianNetwork([ContingencyTable(DataFrame(quote i = [1,2]
    probab = [0.4,0.6] end)), ContingencyTable(DataFrame(quote i = [1,1,2,2]
    j = [1,2,1,2]
    probab = [0.2,0.8,0.3,0.7] end),Array((Function, String),0),["i"])])
    
    wiki_ex = BayesianNetwork([ContingencyTable(DataFrame(quote R = [1,0]
    probab = [0.2,0.8] end)),ContingencyTable(DataFrame(quote R = [0,0,1,1]
    S = [0,1,0,1]
    probab = [0.6,0.4,0.99,0.01] end),Array((Function, String),0),["R"]),ContingencyTable(DataFrame(quote R = [0,0,0,0,1,1,1,1]
    S = [0,0,1,1,0,0,1,1]
    G = [0,1,0,1,0,1,0,1]
    probab = [1,0,0.1,0.9,0.2,0.8,0.01,0.99] end),Array((Function, String),0),["R","S"])])
    
    lecture_ex = BayesianNetwork([ContingencyTable(DataFrame(quote girlfriend = [1,0]
    probab = [0.3,0.7] end)), ContingencyTable(DataFrame(quote weather = ["G","B"]
    probab = [0.25, 0.75] end)), ContingencyTable(DataFrame(quote prof = ["S","B","T"]
    probab = [0.3,0.2,0.5] end)), ContingencyTable(DataFrame(quote topic = ["ML","ML","ML","DB","DB","DB"]
    prof = ["S","B","T","S","B","T"]
    probab = [0.2,0.6,0.9,0.8,0.4,0.1] end),Array((Function,String),0),["prof"]), ContingencyTable(DataFrame(quote girlfriend = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    weather = ["G","G","G","G","G","G","G","G","G","G","G","G","B","B","B","B","B","B","B","B","B","B","B","B","G","G","G","G","G","G","G","G","G","G","G","G","B","B","B","B","B","B","B","B","B","B","B","B"]
    prof = ["S","S","S","S","B","B","B","B","T","T","T","T","S","S","S","S","B","B","B","B","T","T","T","T","S","S","S","S","B","B","B","B","T","T","T","T","S","S","S","S","B","B","B","B","T","T","T","T"]
    topic = ["ML","ML","DB","DB","ML","ML","DB","DB","ML","ML","DB","DB","ML","ML","DB","DB","ML","ML","DB","DB","ML","ML","DB","DB","ML","ML","DB","DB","ML","ML","DB","DB","ML","ML","DB","DB","ML","ML","DB","DB","ML","ML","DB","DB","ML","ML","DB","DB"]
    attend = [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]
    probab = [0.1,0.9,0.2,0.8,0.1,0.9,0.05,0.95,0.3,0.7,0.15,0.85,0.2,0.8,0.25,0.75,0.15,0.85,0.1,0.9,0.35,0.65,0.2,0.8,0.4,0.6,0.6,0.4,0.3,0.7,0.2,0.8,0.7,0.3,0.55,0.45,0.85,0.15,0.9,0.1,0.6,0.4,0.5,0.5,0.95,0.05,0.85,0.15] end),Array((Function,String),0),["prof","topic","girlfriend","weather"])])
    
    test_ex = BayesianNetwork([ContingencyTable(DataFrame(quote A = [0,1]
    probab = [0.4,0.6] end)),ContingencyTable(DataFrame(quote A = [1,1,0,0]
    B = [1,0,1,0]
    probab = [0.1,0.9,0.7,0.3] end),Array((Function,String),0), ["A"]), ContingencyTable(DataFrame(quote C = [0,1]
    probab = [0.2,0.8] end)), ContingencyTable(DataFrame(quote C = [0,0,0,0,1,1,1,1]
    B = [0,0,1,1,0,0,1,1]
    D = [0,1,0,1,0,1,0,1]
    probab = [0.55,0.45,0.8,0.2,0.3,0.7,0.1,0.9] end), Array((Function,String),0),["B","C"]), ContingencyTable(DataFrame(quote D = [0,0,0,0,1,1,1,1]
    B = [0,0,1,1,0,0,1,1]
    E = [0,1,0,1,0,1,0,1]
    probab = [0.65,0.35,0.25,0.75,0.3,0.7,0.05,0.95] end),Array((Function,String),0),["B","D"]),ContingencyTable(DataFrame(quote E = [0,0,0,0,1,1,1,1]
    B = [0,0,1,1,0,0,1,1]
    F = [0,1,0,1,0,1,0,1]
    probab = [0.2,0.8,0.3,0.7,0.4,0.6,0.5,0.5] end),Array((Function,String),0),["B","E"]), ContingencyTable(DataFrame(quote F = [0,0,1,1]
    G = [0,1,0,1]
    probab = [0.1,0.9,0.8,0.2] end),Array((Function,String),0),["F"]),ContingencyTable(DataFrame(quote E = [0,0,1,1]
    H = [0,1,0,1]
    probab = [0.4,0.6,0.8,0.2] end),Array((Function,String),0),["E"]),ContingencyTable(DataFrame(quote H = [0,0,0,0,1,1,1,1]
    C = [0,0,1,1,0,0,1,1]
    I = [0,1,0,1,0,1,0,1]
    probab = [0.25,0.75,0.15,0.85,0.5,0.5,0.6,0.4] end),Array((Function,String),0),["H","C"])])
    
       
    println(out_neighbors(vertices(a.net)[1],a.net))
    println(varnames(a))
    println(childs(a,vertices(a.net)[1]))
    println(parents(a,"i"))
    println(joinedprob(a))
    println(marginal(a,"i"))
    println(marginal(a,"i",[("j",x->x==2)]))
    println(marginal(wiki_ex,"R",[("G",x->x==1)]))
    println(marginal(wiki_ex,"R"))
    tic()
    println(marginal(lecture_ex, "prof",[("attend",x->x==0),("weather",x->x=="B"),("girlfriend",x->x==0)]))
    toc()
    tic()
    println(varelim(lecture_ex, "prof",[("attend",x->x==0),("weather",x->x=="B"),("girlfriend",x->x==0)]))
    toc()
    #println(inducedgraph(lecture_ex))
    tic()
    println(marginal(test_ex, "A",[("D",x->x==0)]))
    toc()
    tic()
    println(varelim(test_ex, "A",[("D",x->x==0)]))
    toc()
    #println(marginal(small_lecture_ex,"topic"))
end