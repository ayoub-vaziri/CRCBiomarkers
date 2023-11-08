# CRCTest

CRCTest is a diagnostic tool based on gene biomarkers for colorectal cancer diagnosis

## Installation

```R
if (!require("randomForest", quietly = T)) install.packages("randomForest")
```

## Usage

```bash
cd ~/
mkdir CRC
cd CRC/
cp YOUR_PATH/CRCTest/RF.rds .
cp YOUR_PATH/CRCTest/RF.R .
cp YOUR_PATH/CRCTest/CRCTest .
chmod +x CRCTest
./CRCTest
```
