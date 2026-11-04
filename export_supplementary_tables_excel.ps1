param(
    [string]$InputDirectory = (Join-Path $PSScriptRoot "supplementary_tables"),
    [string]$OutputDirectory = (Join-Path $PSScriptRoot "supplementary_tables\excel")
)

$ErrorActionPreference = "Stop"

$xlOpenXMLWorkbook = 51
$xlContinuous = 1
$xlThin = 2
$xlEdgeTop = 8
$xlEdgeBottom = 9
$xlCenter = -4108
$xlLeft = -4131
$xlLandscape = 2
$xlPaperA4 = 9

function Convert-CellValue {
    param([AllowNull()][string]$Value)

    if ([string]::IsNullOrWhiteSpace($Value)) {
        return $null
    }

    if ($Value -eq "TRUE") {
        return $true
    }
    if ($Value -eq "FALSE") {
        return $false
    }

    $number = 0.0
    $styles = [System.Globalization.NumberStyles]::Float
    $culture = [System.Globalization.CultureInfo]::InvariantCulture
    if ([double]::TryParse($Value, $styles, $culture, [ref]$number)) {
        return $number
    }

    # Prevent Excel from converting sensitivity counts such as "7/7" to dates.
    if ($Value -match "^\d+\/\d+$") {
        return "'$Value"
    }

    return $Value
}

function Release-ComObject {
    param([AllowNull()]$Object)

    if ($null -ne $Object -and [System.Runtime.InteropServices.Marshal]::IsComObject($Object)) {
        [void][System.Runtime.InteropServices.Marshal]::FinalReleaseComObject($Object)
    }
}

function Export-CompactS10Workbook {
    param(
        [Parameter(Mandatory = $true)]$Excel,
        [Parameter(Mandatory = $true)][array]$Rows,
        [Parameter(Mandatory = $true)][string]$OutputPath
    )

    # Table S10 combines diagnostics that are meaningful for different model
    # structures. Stacking them avoids a wide matrix dominated by blank cells.
    $columnCount = 7
    $rowCount = 34
    $values = [object[,]]::new($rowCount, $columnCount)

    $values[0, 0] = "Table S10. Diagnostic summary for hierarchical cyclic GAMs"
    $values[2, 0] = "Panel A. Sample and model structure"
    $panelAHeaders = @("Model", "Rows", "Camera-rounds", "Trials", "Positive camera-day-hours", "k", "Hierarchy")
    for ($column = 0; $column -lt $panelAHeaders.Count; $column++) {
        $values[3, $column] = $panelAHeaders[$column]
    }
    for ($row = 0; $row -lt $Rows.Count; $row++) {
        $values[(4 + $row), 0] = $Rows[$row].model
        $values[(4 + $row), 1] = Convert-CellValue ([string]$Rows[$row].aggregated_rows)
        $values[(4 + $row), 2] = Convert-CellValue ([string]$Rows[$row].camera_rounds)
        $values[(4 + $row), 3] = Convert-CellValue ([string]$Rows[$row].binomial_trials)
        $values[(4 + $row), 4] = Convert-CellValue ([string]$Rows[$row].positive_camera_day_hours)
        $values[(4 + $row), 5] = Convert-CellValue ([string]$Rows[$row].cyclic_basis_k)
        $values[(4 + $row), 6] = $Rows[$row].hierarchy
    }

    $values[10, 0] = "Panel B. Common fitting diagnostics"
    $panelBHeaders = @("Model", "Deviance explained", "Pearson dispersion", "Convergence")
    for ($column = 0; $column -lt $panelBHeaders.Count; $column++) {
        $values[11, $column] = $panelBHeaders[$column]
    }
    for ($row = 0; $row -lt $Rows.Count; $row++) {
        $values[(12 + $row), 0] = $Rows[$row].model
        $values[(12 + $row), 1] = Convert-CellValue ([string]$Rows[$row].deviance_explained)
        $values[(12 + $row), 2] = Convert-CellValue ([string]$Rows[$row].Pearson_dispersion)
        $values[(12 + $row), 3] = $Rows[$row].convergence
    }

    $values[18, 0] = "Panel C. Structure-specific diagnostic checks"
    $panelCHeaders = @("Model", "Diagnostic", "Value")
    for ($column = 0; $column -lt $panelCHeaders.Count; $column++) {
        $values[19, $column] = $panelCHeaders[$column]
    }

    $classRows = @($Rows | Where-Object { $_.model -ne "Confirmatory joint GAM" })
    $juvenile = @($Rows | Where-Object { $_.model -eq "Juvenile" })[0]
    $joint = @($Rows | Where-Object { $_.model -eq "Confirmatory joint GAM" })[0]
    $specificChecks = New-Object System.Collections.ArrayList
    foreach ($row in $classRows) {
        [void]$specificChecks.Add([object[]]@($row.model, "Mean camera-round lag-1 Pearson residual correlation", $row.mean_camera_round_lag1_Pearson))
    }
    [void]$specificChecks.Add([object[]]@("Juvenile", "Minimum time-smooth k-index", $juvenile.minimum_time_k_index))
    [void]$specificChecks.Add([object[]]@("Juvenile", "Time-smooth k-test P", $juvenile.minimum_time_k_P))
    [void]$specificChecks.Add([object[]]@("Juvenile", "Maximum time-smooth concurvity", $juvenile.maximum_time_concurvity))
    [void]$specificChecks.Add([object[]]@("Confirmatory joint GAM", "Maximum absolute gradient", $joint.maximum_absolute_gradient))
    [void]$specificChecks.Add([object[]]@("Confirmatory joint GAM", "Minimum Hessian eigenvalue", $joint.minimum_hessian_eigenvalue))
    [void]$specificChecks.Add([object[]]@("Confirmatory joint GAM", "Hessian numerically positive semidefinite", $joint.hessian_numerically_PSD))
    [void]$specificChecks.Add([object[]]@("Confirmatory joint GAM", "Coefficient rank", $joint.coefficient_rank))

    for ($row = 0; $row -lt $specificChecks.Count; $row++) {
        $values[(20 + $row), 0] = $specificChecks[$row][0]
        $values[(20 + $row), 1] = $specificChecks[$row][1]
        $values[(20 + $row), 2] = Convert-CellValue ([string]$specificChecks[$row][2])
    }

    $values[32, 0] = "Structure-specific diagnostics are listed only where reported or meaningful. The juvenile k-index motivated the prespecified k = 8 and k = 16 sensitivity fits; no explicit AR(1) correction was used."
    $values[33, 0] = "Sources: definitive_hierarchical_diel_diagnostics.csv; definitive_hierarchical_diel_data_audit.csv; joint_demographic_penalty_calibrated_outputs/penalty_calibrated_diagnostics.csv."

    $workbook = $null
    $worksheet = $null
    $usedRange = $null
    try {
        $workbook = $Excel.Workbooks.Add()
        $worksheet = $workbook.Worksheets.Item(1)
        $worksheet.Name = "Table S10"

        $usedRange = $worksheet.Range($worksheet.Cells.Item(1, 1), $worksheet.Cells.Item($rowCount, $columnCount))
        $usedRange.Value2 = $values
        $usedRange.Font.Name = "Times New Roman"
        $usedRange.Font.Size = 10
        $usedRange.VerticalAlignment = $xlCenter
        $usedRange.WrapText = $true

        foreach ($mergeRow in @(1, 3, 11, 19, 33, 34)) {
            $mergeRange = $worksheet.Range($worksheet.Cells.Item($mergeRow, 1), $worksheet.Cells.Item($mergeRow, $columnCount))
            [void]$mergeRange.Merge()
            Release-ComObject $mergeRange
        }

        $titleRange = $worksheet.Range($worksheet.Cells.Item(1, 1), $worksheet.Cells.Item(1, $columnCount))
        $titleRange.Font.Bold = $true
        $titleRange.Font.Size = 11
        $titleRange.HorizontalAlignment = $xlLeft
        Release-ComObject $titleRange

        foreach ($panelRow in @(3, 11, 19)) {
            $panelRange = $worksheet.Range($worksheet.Cells.Item($panelRow, 1), $worksheet.Cells.Item($panelRow, $columnCount))
            $panelRange.Font.Bold = $true
            $panelRange.HorizontalAlignment = $xlLeft
            Release-ComObject $panelRange
        }

        foreach ($headerRow in @(4, 12, 20)) {
            $headerRange = $worksheet.Range($worksheet.Cells.Item($headerRow, 1), $worksheet.Cells.Item($headerRow, $columnCount))
            $headerRange.Font.Bold = $true
            $headerRange.HorizontalAlignment = $xlCenter
            foreach ($edge in @($xlEdgeTop, $xlEdgeBottom)) {
                $border = $headerRange.Borders.Item($edge)
                $border.LineStyle = $xlContinuous
                $border.Weight = $xlThin
                $border.Color = 0
                Release-ComObject $border
            }
            Release-ComObject $headerRange
        }

        foreach ($lastDataRow in @(9, 17, 31)) {
            $lastRowRange = $worksheet.Range($worksheet.Cells.Item($lastDataRow, 1), $worksheet.Cells.Item($lastDataRow, $columnCount))
            $border = $lastRowRange.Borders.Item($xlEdgeBottom)
            $border.LineStyle = $xlContinuous
            $border.Weight = $xlThin
            $border.Color = 0
            Release-ComObject $border
            Release-ComObject $lastRowRange
        }

        $worksheet.Columns.Item(1).ColumnWidth = 24
        $worksheet.Columns.Item(2).ColumnWidth = 31
        $worksheet.Columns.Item(3).ColumnWidth = 17
        $worksheet.Columns.Item(4).ColumnWidth = 18
        $worksheet.Columns.Item(5).ColumnWidth = 22
        $worksheet.Columns.Item(6).ColumnWidth = 8
        $worksheet.Columns.Item(7).ColumnWidth = 46
        [void]$usedRange.Rows.AutoFit()
        $worksheet.Rows.Item(1).RowHeight = 22

        $footnoteRange = $worksheet.Range($worksheet.Cells.Item(33, 1), $worksheet.Cells.Item(34, $columnCount))
        $footnoteRange.Font.Size = 9
        $footnoteRange.Font.Italic = $true
        $footnoteRange.HorizontalAlignment = $xlLeft
        Release-ComObject $footnoteRange

        [void]$worksheet.Activate()
        $Excel.ActiveWindow.DisplayGridlines = $false
        $worksheet.PageSetup.Orientation = $xlLandscape
        $worksheet.PageSetup.PaperSize = $xlPaperA4
        $worksheet.PageSetup.Zoom = $false
        $worksheet.PageSetup.FitToPagesWide = 1
        $worksheet.PageSetup.FitToPagesTall = $false
        $worksheet.PageSetup.CenterHorizontally = $true
        $worksheet.PageSetup.LeftMargin = $Excel.InchesToPoints(0.4)
        $worksheet.PageSetup.RightMargin = $Excel.InchesToPoints(0.4)
        $worksheet.PageSetup.TopMargin = $Excel.InchesToPoints(0.5)
        $worksheet.PageSetup.BottomMargin = $Excel.InchesToPoints(0.5)
        $worksheet.PageSetup.PrintArea = $usedRange.Address()

        $workbook.SaveAs($OutputPath, $xlOpenXMLWorkbook)
    }
    finally {
        if ($null -ne $workbook) {
            $workbook.Close($false)
        }
        Release-ComObject $usedRange
        Release-ComObject $worksheet
        Release-ComObject $workbook
        [GC]::Collect()
        [GC]::WaitForPendingFinalizers()
    }
}

$inputRoot = (Resolve-Path -LiteralPath $InputDirectory).Path
if (-not (Test-Path -LiteralPath $OutputDirectory)) {
    New-Item -ItemType Directory -Path $OutputDirectory | Out-Null
}
$outputRoot = (Resolve-Path -LiteralPath $OutputDirectory).Path

$files = @(Get-ChildItem -LiteralPath $inputRoot -Filter "table_S*.csv" -File | Sort-Object Name)
if ($files.Count -ne 10) {
    throw "Expected exactly 10 supplementary CSV tables; found $($files.Count)."
}

$excel = $null
try {
    $excel = New-Object -ComObject Excel.Application
    $excel.Visible = $false
    $excel.DisplayAlerts = $false
    $excel.ScreenUpdating = $false
    $excel.SheetsInNewWorkbook = 1

    foreach ($file in $files) {
        $rows = @(Import-Csv -LiteralPath $file.FullName)
        if ($rows.Count -eq 0) {
            throw "Table has no data rows: $($file.FullName)"
        }

        $tableNumber = [int]([regex]::Match($file.BaseName, "S(\d+)").Groups[1].Value)
        $outputName = [System.IO.Path]::ChangeExtension($file.Name, ".xlsx")
        $outputPath = [System.IO.Path]::GetFullPath((Join-Path $outputRoot $outputName))
        if (-not $outputPath.StartsWith($outputRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
            throw "Refusing to write outside the intended output directory: $outputPath"
        }

        if ($tableNumber -eq 10) {
            Export-CompactS10Workbook -Excel $excel -Rows $rows -OutputPath $outputPath
            Write-Output $outputPath
            continue
        }

        $headers = @($rows[0].PSObject.Properties.Name)
        $rowCount = $rows.Count + 1
        $columnCount = $headers.Count
        $values = [object[,]]::new($rowCount, $columnCount)

        for ($column = 0; $column -lt $columnCount; $column++) {
            $values[0, $column] = $headers[$column]
        }
        for ($row = 0; $row -lt $rows.Count; $row++) {
            for ($column = 0; $column -lt $columnCount; $column++) {
                $raw = [string]$rows[$row].PSObject.Properties[$headers[$column]].Value
                $values[($row + 1), $column] = Convert-CellValue $raw
            }
        }

        $workbook = $null
        $worksheet = $null
        $usedRange = $null
        $headerRange = $null
        $lastRowRange = $null
        try {
            $workbook = $excel.Workbooks.Add()
            $worksheet = $workbook.Worksheets.Item(1)

            $worksheet.Name = "Table S$tableNumber"

            $topLeft = $worksheet.Cells.Item(1, 1)
            $bottomRight = $worksheet.Cells.Item($rowCount, $columnCount)
            $usedRange = $worksheet.Range($topLeft, $bottomRight)
            $usedRange.Value2 = $values
            $usedRange.Font.Name = "Times New Roman"
            $usedRange.Font.Size = 10
            $usedRange.VerticalAlignment = $xlCenter

            $headerRange = $worksheet.Range(
                $worksheet.Cells.Item(1, 1),
                $worksheet.Cells.Item(1, $columnCount)
            )
            $headerRange.Font.Bold = $true
            $headerRange.Font.Size = 10
            $headerRange.HorizontalAlignment = $xlCenter
            $headerRange.WrapText = $true

            foreach ($edge in @($xlEdgeTop, $xlEdgeBottom)) {
                $border = $headerRange.Borders.Item($edge)
                $border.LineStyle = $xlContinuous
                $border.Weight = $xlThin
                $border.Color = 0
                Release-ComObject $border
            }

            $lastRowRange = $worksheet.Range(
                $worksheet.Cells.Item($rowCount, 1),
                $worksheet.Cells.Item($rowCount, $columnCount)
            )
            $bottomBorder = $lastRowRange.Borders.Item($xlEdgeBottom)
            $bottomBorder.LineStyle = $xlContinuous
            $bottomBorder.Weight = $xlThin
            $bottomBorder.Color = 0
            Release-ComObject $bottomBorder

            $dataRange = $worksheet.Range(
                $worksheet.Cells.Item(2, 1),
                $worksheet.Cells.Item($rowCount, $columnCount)
            )
            $dataRange.HorizontalAlignment = $xlLeft
            $dataRange.WrapText = $true

            for ($column = 1; $column -le $columnCount; $column++) {
                $columnRange = $worksheet.Columns.Item($column)
                [void]$columnRange.AutoFit()

                $headerName = $headers[$column - 1]
                $maximumWidth = if ($headerName -match "source|note|formula|structure|hierarchy|interpretation|contrast") {
                    42
                } else {
                    24
                }
                if ($columnRange.ColumnWidth -gt $maximumWidth) {
                    $columnRange.ColumnWidth = $maximumWidth
                }
                if ($columnRange.ColumnWidth -lt 8) {
                    $columnRange.ColumnWidth = 8
                }
                Release-ComObject $columnRange
            }

            [void]$usedRange.Rows.AutoFit()
            if ($worksheet.Rows.Item(1).RowHeight -lt 30) {
                $worksheet.Rows.Item(1).RowHeight = 30
            }

            [void]$worksheet.Activate()
            $excel.ActiveWindow.DisplayGridlines = $false
            $excel.ActiveWindow.SplitRow = 1
            $excel.ActiveWindow.FreezePanes = $true

            $worksheet.PageSetup.Orientation = $xlLandscape
            $worksheet.PageSetup.PaperSize = $xlPaperA4
            $worksheet.PageSetup.Zoom = $false
            $worksheet.PageSetup.FitToPagesWide = 1
            $worksheet.PageSetup.FitToPagesTall = $false
            $worksheet.PageSetup.CenterHorizontally = $true
            $worksheet.PageSetup.LeftMargin = $excel.InchesToPoints(0.4)
            $worksheet.PageSetup.RightMargin = $excel.InchesToPoints(0.4)
            $worksheet.PageSetup.TopMargin = $excel.InchesToPoints(0.5)
            $worksheet.PageSetup.BottomMargin = $excel.InchesToPoints(0.5)
            $worksheet.PageSetup.PrintTitleRows = '$1:$1'
            $worksheet.PageSetup.PrintArea = $usedRange.Address()

            $workbook.SaveAs($outputPath, $xlOpenXMLWorkbook)
            Write-Output $outputPath
        }
        finally {
            if ($null -ne $workbook) {
                $workbook.Close($false)
            }
            Release-ComObject $lastRowRange
            Release-ComObject $headerRange
            Release-ComObject $usedRange
            Release-ComObject $worksheet
            Release-ComObject $workbook
            [GC]::Collect()
            [GC]::WaitForPendingFinalizers()
        }
    }
}
finally {
    if ($null -ne $excel) {
        $excel.Quit()
    }
    Release-ComObject $excel
    [GC]::Collect()
    [GC]::WaitForPendingFinalizers()
}
