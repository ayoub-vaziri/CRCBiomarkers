# CRCTest

**CRCTest is a diagnostic tool based on gene biomarkers for colorectal cancer diagnosis**

## Installation

```R
if (!require("randomForest", quietly = T)) 
	install.packages("randomForest")
```

- require("randomForest", quietly = TRUE): This checks if the "randomForest" package is already installed. If it is, it returns TRUE; otherwise, it returns FALSE.
- !: The exclamation mark is used to negate the result. So, !TRUE becomes FALSE, and !FALSE becomes TRUE.
- install.packages("randomForest"): If the package is not installed (require returns FALSE), this line installs the "randomForest" package.

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