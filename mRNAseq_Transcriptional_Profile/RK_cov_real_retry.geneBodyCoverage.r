accepted_hits <- c(0.0,0.09535673839184598,0.16187869636340757,0.20387488989555808,0.24689977349943376,0.2756912986032465,0.3016720145967032,0.33097945765697745,0.35296338240845604,0.3776487982886624,0.4132597206493016,0.4297478608279854,0.4517168428337737,0.46217991065811,0.4907779665282497,0.5153320435384422,0.527320057883478,0.5360356109223606,0.5597277274443186,0.5811697810494526,0.6010585755631056,0.617464923870643,0.639256165848748,0.6501242607273184,0.6687429218573047,0.6835047816786208,0.6941573864351328,0.7203142695356738,0.72388715867623,0.7235521265886498,0.7136678935447338,0.7269881716370957,0.7351657858311312,0.7533676230023908,0.7692407512268781,0.7769173902101422,0.78233216937209,0.8129899647665786,0.8080140304517428,0.814524348810872,0.8017836919592299,0.8092275386938468,0.8206713225116397,0.830216591166478,0.8466339499182082,0.8479866616333207,0.8510286900717252,0.8758234239335598,0.892193595067321,0.9035524411727696,0.9001226878067196,0.8991694979237448,0.9007550018875047,0.9162569208506355,0.9126014533786334,0.9128428966905751,0.9061886560966402,0.9082814269535674,0.9040479111614446,0.9209017553793885,0.9232564175160438,0.9288961243236441,0.9429855605889015,0.9423233610167359,0.9652596891908897,0.9799995281238203,0.9842723669309174,0.9837690323392475,0.969614319869133,0.9569814080785202,0.9595515603372342,0.9493621806971184,0.9632361268403171,0.9615782685290046,0.95207310934944,0.9603427393985151,0.9693083868126336,0.973144740153517,0.9888590033975085,1.0,0.9944987102051088,0.9816597458160312,0.9727444318610796,0.9703205612180698,0.964750062916824,0.9649985843714609,0.9423406631433245,0.9374087706052598,0.9335071410595193,0.9115601799421166,0.8935423744809362,0.8952670819177048,0.8791989115389455,0.8619652069963508,0.825248521454637,0.7880410846860451,0.7459127658235812,0.6996232855165472,0.6143969422423556,0.4353883540958852)


pdf("RK_cov_real_retry.geneBodyCoverage.curves.pdf")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot(x,accepted_hits,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
dev.off()
