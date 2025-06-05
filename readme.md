# 🇳🇱 NL 2030 V2G Analysis — Merit Order & Battery Scenario Simulation

This MATLAB script performs a simulation of the Dutch electricity market in the year 2030 under different **Vehicle-to-Grid (V2G)** battery storage scenarios. Developed as part of the master's thesis by **Mayk Thewessen**, the model explores how increasing V2G capacity affects the **hourly electricity price profile**, based on **merit order dispatch**, **solar PV** and **wind generation assumptions**, and projected **demand profiles**.

## 🧠 Objective

To assess the impact of V2G penetration in the Netherlands in 2030 on:
- Hourly electricity prices
- Load shifting potential
- Renewable curtailment avoidance
- System-wide efficiency based on dynamic dispatch

## 🧮 Core Concepts

- **Merit Order Pricing:** Simulated hourly market prices based on generator marginal cost ranking.
- **V2G Batteries:** Simulated as flexible demand/supply depending on grid condition.
- **Scenario Years:** Simulation can span multiple test years.
- **Renewables in 2030:** Includes predefined solar PV and wind capacity.

## 📂 Main Steps in the Script

1. **Initialize Environment** (clear memory, set format)
2. **Load Demand Data** (`readtable` from XML/XLSX/CSV formats)
3. **Reshape & Convert Load Vectors**
4. **Define V2G Storage & Renewable Scenarios**
5. **Loop Over Simulation Years**
6. **Calculate Hourly Merit Order**
7. **Adjust Load/Price with V2G Dispatch**
8. **Export Results to Excel**

## 📊 Example Output

- Excel files with hourly price profiles
- Time-series plots of load vs. generation
- Scenario comparison of price volatility

## 📈 Screenshot

![V2G Impact Simulation](./preview.svg)

## 📦 Dependencies

- MATLAB R2020b or higher recommended
- Compatible with `.xlsx`, `.xml`, and `.csv` data input
- Uses built-in `readtable`, `plot`, and `xlswrite`

## 📘 File Overview

| File                        | Description                                      |
|-----------------------------|--------------------------------------------------|
| `Lipton_v4_6_export_xlsx.m` | Main simulation script                          |
| `ACTUAL_TOTAL_LOAD_*.xml`   | Input files with hourly load data               |
| `*.xlsx` or `*.csv`         | (Optional) Alternative data input formats       |

## 🧪 Example Use

```matlab
run('Lipton_v4_6_export_xlsx.m');