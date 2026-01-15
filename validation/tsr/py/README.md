# SHUD .dat Reader (Python)

This folder contains a small Python reader for SHUD binary output files (`.dat`).

## File format (SHUD `.dat`)

Based on `src/classes/Model_Control.cpp` (`Print_Ctrl::open_file`), the layout is:

1. **1024 bytes**: fixed-size ASCII header (padded with `\\0`)
2. **8 bytes**: `StartTime` (C++ `double`)
3. **8 bytes**: `NumVar` (C++ `double`, integer-valued)
4. **8 * NumVar bytes**: `icol` array (C++ `double[]`)
5. **Payload**: `NumRecords` repeated records, each record is:
   - `Time_min` (double)
   - `NumVar` values (double[])

## Install deps

`DatReader.to_dataframe()` requires `pandas` (and `numpy`, via pandas):

```bash
python3 -m pip install pandas
```

## API

- `validation/tsr/py/shud_reader.py`
  - `DatReader(path)`: parses header + reads payload
  - `DatReader.meta`: parsed metadata (`DatMeta`)
  - `DatReader.to_dataframe(time_mode="index"|"column")`
  - `read_dat(path, ...)`: convenience wrapper returning a DataFrame

Example:

```python
from shud_reader import DatReader

r = DatReader("output/ccw.base/ccw.rn_factor.dat")
df = r.to_dataframe(time_mode="index")  # index is Time_min
print(df.shape)
```

## CLI helper

List `.dat` files and basic metadata:

```bash
python3 validation/tsr/py/inspect_output.py output/ccw.base
```

Export one `.dat` to CSV:

```bash
python3 validation/tsr/py/inspect_output.py output/ccw.base --export ccw.rn_factor.dat
```

## Tests

Run unit tests:

```bash
python3 -m unittest discover -s validation/tsr/py -p 'test_*.py'
```

To enforce coverage (requires `coverage`):

```bash
python3 -m pip install coverage
coverage run -m unittest discover -s validation/tsr/py -p 'test_*.py'
coverage report --fail-under=90
```

