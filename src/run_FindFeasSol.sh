m=(1000 2000 4000 6000 8000 10000 20000 40000 60000 80000 100000)
n=(50 100)
for i in "${m[@]}"; do
  for j in "${n[@]}"; do
    echo "Running FindFeasSol.py with argument $i $j"
    python3 FindFeasSol.py "$i" "$j"
    echo "------------------------"
  done  
done
