# Query Parameter API Documentation

A URL-based configuration system for sharing application state via query parameters.

## Overview

This system allows you to configure the application by appending parameters to the URL. All parameters are optional and have sensible defaults.

## Base URL Format

```
https://reverse-metabolomics.gnps2.org/?param1=value1&param2=value2&...
```

## Available Parameters

### `usi` (JSON Object)
**Type:** JSON  
**Default:** None  
**Description:** USI (Universal Spectrum Identifier) data mapping compound names to spectrum identifiers.

**Format:**
```
?usi={"compound_name":"mzspec:...","another_compound":"mzspec:..."}
```

**Example:**
```
?usi=%7B%22Caffeine%22%3A%22mzspec%3AGNPS%3ATASK-123%22%7D
```

**Notes:**
- URL-encoded JSON object
- Keys: compound names
- Values: mzspec identifiers
- Mutually exclusive with `example` parameter

---

### `example` (String)
**Type:** String (enum)  
**Default:** None  
**Valid Values:** `"1"`, `"2"`  
**Description:** Load a predefined example dataset.

**Example:**
```
?example=1
```

**Notes:**
- Setting this parameter removes `usi` from URL
- Use to quickly load demo data

---

### `database` (String)
**Type:** String  
**Default:** `"metabolomicspanrepo_index_latest"`  
**Description:** Specifies which MASST database to search.

**Example:**
```
?database=gnpsdata_index
```

---

### `analog` (Boolean)
**Type:** Boolean  
**Default:** `false`  
**Valid Values:** `true`, `false`, `1`, `0`, `yes`, `no`  
**Description:** Enable analog search mode.

**Example:**
```
?analog=true
```

---

### `precursor_tol` (Float)
**Type:** Float  
**Default:** `0.02`  
**Range:** `0.001` to `1.0`  
**Description:** Precursor m/z tolerance in Daltons.

**Example:**
```
?precursor_tol=0.05
```

---

### `fragment_tol` (Float)
**Type:** Float  
**Default:** `0.02`  
**Range:** `0.001` to `1.0`  
**Description:** Fragment m/z tolerance in Daltons.

**Example:**
```
?fragment_tol=0.01
```

---

### `min_cos` (Float)
**Type:** Float  
**Default:** `0.7`  
**Range:** `0.0` to `1.0`  
**Description:** Minimum cosine similarity threshold for matches.

**Example:**
```
?min_cos=0.8
```

---

### `organism` (String)
**Type:** String (enum)  
**Default:** `"Humans"`  
**Valid Values:** `"Humans"`, `"Rodents"`, `"All organisms"`, `"Manual Entry"`  
**Description:** Filter results by organism type.

**Example:**
```
?organism=Rodents
```

---

## Complete Example URLs

### Example 1: Basic search with custom tolerances
```
https://reverse-metabolomics.gnps2.org/?database=gnpsdata_index&precursor_tol=0.05&fragment_tol=0.03&min_cos=0.75
```

### Example 2: Analog search with specific organism
```
https://reverse-metabolomics.gnps2.org/?analog=true&organism=Rodents&min_cos=0.6
```

### Example 3: Load example dataset
```
https://reverse-metabolomics.gnps2.org/?example=1&database=metabolomicspanrepo_index_latest
```

### Example 4: Custom USI data
```
https://reverse-metabolomics.gnps2.org/?usi=%7B%22Caffeine%22%3A%22mzspec%3AGNPS%3ATASK-4f69e11bfb544010b2c4117fac34a5b5-spectra%2Fspecs_ms.mgf%3Ascan%3A1840%22%7D&min_cos=0.8
```

## Parameter Schema Export

To programmatically discover available parameters, call:

```python
from query_params import get_param_schema

schema = get_param_schema()
# Returns list of OpenAPI-style parameter definitions
```

## API Functions

### Loading Parameters
```python
from query_params import load_all_params

params = load_all_params()
# Returns: {'usi': None, 'example': None, 'database': 'metabolomicspanrepo_index_latest', ...}
```

### Syncing Parameters
```python
from query_params import sync_params, sync_param

# Sync multiple parameters
sync_params({
    'analog': True,
    'min_cos': 0.8,
    'organism': 'Rodents'
})

# Sync single parameter
sync_param('precursor_tol', 0.05)
```

### Getting Defaults
```python
from query_params import get_defaults, get_default

# Get all defaults
defaults = get_defaults()

# Get specific default
default_tol = get_default('precursor_tol')  # Returns: 0.02
```

## Notes

- Parameters equal to their default values are omitted from the URL
- Invalid values fall back to defaults
- All numeric parameters are validated against their allowed ranges
- Boolean parameters accept various formats: `true/false`, `1/0`, `yes/no`
- JSON parameters are URL-encoded automatically
