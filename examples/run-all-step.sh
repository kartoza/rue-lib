#!/bin/bash

# Function to format elapsed time
format_time() {
    local seconds=$1
    printf "%02d:%02d:%02d" $((seconds/3600)) $((seconds%3600/60)) $((seconds%60))
}

# Start total timer
TOTAL_START=$SECONDS

echo "========================================="
echo "Starting all steps"
echo "========================================="

# Step 1
echo ""
echo "Running Step 1: Generate Parcels..."
STEP_START=$SECONDS
python examples/step1_generate_parcels.py
STEP_ELAPSED=$((SECONDS - STEP_START))
echo "Step 1 completed in $(format_time $STEP_ELAPSED)"

# Step 2
echo ""
echo "Running Step 2: Generate Streets..."
STEP_START=$SECONDS
python examples/step2_generate_streets.py
STEP_ELAPSED=$((SECONDS - STEP_START))
echo "Step 2 completed in $(format_time $STEP_ELAPSED)"

# Step 3
echo ""
echo "Running Step 3: Generate Cluster..."
STEP_START=$SECONDS
python examples/step3_generate_cluster.py
STEP_ELAPSED=$((SECONDS - STEP_START))
echo "Step 3 completed in $(format_time $STEP_ELAPSED)"

# Step 4
echo ""
echo "Running Step 4: Generate Public..."
STEP_START=$SECONDS
python examples/step4_generate_public.py
STEP_ELAPSED=$((SECONDS - STEP_START))
echo "Step 4 completed in $(format_time $STEP_ELAPSED)"

# Total time
TOTAL_ELAPSED=$((SECONDS - TOTAL_START))
echo ""
echo "========================================="
echo "All steps completed!"
echo "Total time: $(format_time $TOTAL_ELAPSED)"
echo "========================================="
