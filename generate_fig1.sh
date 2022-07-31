numberOfProcessors=1
numberOfSimulations=1
rows=25
columns=25
neighborhoodThreshold=1.45
recordRate=10
imageReturn=false
imageRecordRate=1 
r=1.0
d=.4
a=.2
g=.11 
y=.75
dt=.1
tf=1 
blobValue=0
gridOption0=0
gridOption2=2

# create directories
grazingFolder=$(python -c "print(str(round($g * 100)))")
thresholdFolder=$(python -c "print(int($neighborhoodThreshold*100))")
mkdir -p 'output'/$rows'x'$columns/'grid'$gridOption0/'grazing'$grazingFolder/'threshold'$thresholdFolder/'images'
mkdir -p 'output'/$rows'x'$columns/'grid'$gridOption2/'grazing'$grazingFolder/'threshold'$thresholdFolder/'images'

# simulate different initial profiles

coralPercent=33
macroalgaePercent=33
gridOption=0
python coralModel_functions.py $numberOfProcessors $numberOfSimulations $coralPercent $macroalgaePercent $gridOption $rows $columns $neighborhoodThreshold $recordRate $imageReturn $imageRecordRate $r $d $a $g $y $dt $tf $blobValue
 
coralPercent=50
macroalgaePercent=25
gridOption=0
python coralModel_functions.py $numberOfProcessors $numberOfSimulations $coralPercent $macroalgaePercent $gridOption $rows $columns $neighborhoodThreshold $recordRate $imageReturn $imageRecordRate $r $d $a $g $y $dt $tf $blobValue
     
coralPercent=33
macroalgaePercent=33
gridOption=2
blobValue=0
python coralModel_functions.py $numberOfProcessors $numberOfSimulations $coralPercent $macroalgaePercent $gridOption $rows $columns $neighborhoodThreshold $recordRate $imageReturn $imageRecordRate $r $d $a $g $y $dt $tf $blobValue

coralPercent=50
macroalgaePercent=25
gridOption=2
blobValue=0
python coralModel_functions.py $numberOfProcessors $numberOfSimulations $coralPercent $macroalgaePercent $gridOption $rows $columns $neighborhoodThreshold $recordRate $imageReturn $imageRecordRate $r $d $a $g $y $dt $tf $blobValue       

# run example simulations

numberOfProcessors=100
numberOfSimulations=100
g=.58
# create directories
grazingFolder=$(python -c "print(str(round($g * 100)))")
thresholdFolder=$(python -c "print(int($neighborhoodThreshold*100))")
mkdir -p 'output'/$rows'x'$columns/'grid'$gridOption0/'grazing'$grazingFolder/'threshold'$thresholdFolder/'images'
# simulate
coralPercent=33
macroalgaePercent=33
tf=101
gridOption=0
python coralModel_functions.py $numberOfProcessors $numberOfSimulations $coralPercent $macroalgaePercent $gridOption $rows $columns $neighborhoodThreshold $recordRate $imageReturn $imageRecordRate $r $d $a $g $y $dt $tf $blobValue
        
            







