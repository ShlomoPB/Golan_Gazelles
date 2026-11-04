param(
  [string]$SourceMarkdown = "conbio_methods_results_raw_package.md",
  [string]$OutputDocx = "Golan_Gazelles_ConBio_Methods_Results.docx"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$workspace = (Get-Location).Path
$sourcePath = Join-Path $workspace $SourceMarkdown
$outputPath = Join-Path $workspace $OutputDocx

if (-not (Test-Path -LiteralPath $sourcePath)) {
  throw "Source markdown not found: $sourcePath"
}

Add-Type -AssemblyName System.IO.Compression
Add-Type -AssemblyName System.IO.Compression.FileSystem

function Escape-XmlText {
  param([AllowNull()][string]$Text)
  if ($null -eq $Text) { return "" }
  return [System.Security.SecurityElement]::Escape($Text)
}

function Strip-MarkdownInline {
  param([string]$Text)
  $out = $Text
  $out = $out -replace '\[([^\]]+)\]\([^)]+\)', '$1'
  $out = $out -replace '`([^`]+)`', '$1'
  $out = $out -replace '\\\(', '('
  $out = $out -replace '\\\)', ')'
  return $out
}

function New-Paragraph {
  param(
    [string]$Text,
    [ValidateSet("Normal", "Title", "Heading1", "Heading2", "Heading3", "Caption", "Bullet")]
    [string]$Style = "Normal"
  )

  $clean = Escape-XmlText (Strip-MarkdownInline $Text)
  $pPr = ""
  $rPr = ""
  $prefix = ""

  switch ($Style) {
    "Title" {
      $pPr = '<w:pPr><w:spacing w:after="240"/><w:jc w:val="center"/></w:pPr>'
      $rPr = '<w:rPr><w:b/><w:sz w:val="32"/></w:rPr>'
    }
    "Heading1" {
      $pPr = '<w:pPr><w:spacing w:before="300" w:after="120"/></w:pPr>'
      $rPr = '<w:rPr><w:b/><w:sz w:val="28"/></w:rPr>'
    }
    "Heading2" {
      $pPr = '<w:pPr><w:spacing w:before="220" w:after="100"/></w:pPr>'
      $rPr = '<w:rPr><w:b/><w:sz w:val="24"/></w:rPr>'
    }
    "Heading3" {
      $pPr = '<w:pPr><w:spacing w:before="180" w:after="80"/></w:pPr>'
      $rPr = '<w:rPr><w:b/><w:i/><w:sz w:val="22"/></w:rPr>'
    }
    "Caption" {
      $pPr = '<w:pPr><w:spacing w:before="80" w:after="180"/></w:pPr>'
      $rPr = '<w:rPr><w:i/><w:sz w:val="20"/></w:rPr>'
    }
    "Bullet" {
      $pPr = '<w:pPr><w:ind w:left="360" w:hanging="180"/></w:pPr>'
      $prefix = "- "
    }
    default {
      $pPr = '<w:pPr><w:spacing w:after="120"/><w:jc w:val="both"/></w:pPr>'
      $rPr = '<w:rPr><w:sz w:val="22"/></w:rPr>'
    }
  }

  $textValue = Escape-XmlText ($prefix + (Strip-MarkdownInline $Text))
  return "<w:p>$pPr<w:r>$rPr<w:t xml:space=`"preserve`">$textValue</w:t></w:r></w:p>"
}

function New-PageBreak {
  return '<w:p><w:r><w:br w:type="page"/></w:r></w:p>'
}

function New-ImageParagraph {
  param(
    [string]$RelId,
    [string]$Name,
    [double]$WidthInches,
    [double]$HeightInches
  )

  $cx = [math]::Round($WidthInches * 914400)
  $cy = [math]::Round($HeightInches * 914400)
  $docPrId = [math]::Abs($Name.GetHashCode())
  if ($docPrId -eq 0) { $docPrId = 1 }
  $escapedName = Escape-XmlText $Name

  return @"
<w:p>
  <w:pPr><w:jc w:val="center"/><w:spacing w:before="120" w:after="80"/></w:pPr>
  <w:r>
    <w:drawing>
      <wp:inline distT="0" distB="0" distL="0" distR="0" xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing">
        <wp:extent cx="$cx" cy="$cy"/>
        <wp:docPr id="$docPrId" name="$escapedName"/>
        <wp:cNvGraphicFramePr>
          <a:graphicFrameLocks noChangeAspect="1" xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"/>
        </wp:cNvGraphicFramePr>
        <a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">
          <a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">
            <pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">
              <pic:nvPicPr>
                <pic:cNvPr id="0" name="$escapedName"/>
                <pic:cNvPicPr/>
              </pic:nvPicPr>
              <pic:blipFill>
                <a:blip r:embed="$RelId" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"/>
                <a:stretch><a:fillRect/></a:stretch>
              </pic:blipFill>
              <pic:spPr>
                <a:xfrm>
                  <a:off x="0" y="0"/>
                  <a:ext cx="$cx" cy="$cy"/>
                </a:xfrm>
                <a:prstGeom prst="rect"><a:avLst/></a:prstGeom>
              </pic:spPr>
            </pic:pic>
          </a:graphicData>
        </a:graphic>
      </wp:inline>
    </w:drawing>
  </w:r>
</w:p>
"@
}

function New-WordTable {
  param(
    [string[]]$Headers,
    [object[]]$Rows
  )

  $xml = '<w:tbl><w:tblPr><w:tblW w:w="0" w:type="auto"/><w:tblBorders><w:top w:val="single" w:sz="4" w:space="0" w:color="999999"/><w:left w:val="single" w:sz="4" w:space="0" w:color="999999"/><w:bottom w:val="single" w:sz="4" w:space="0" w:color="999999"/><w:right w:val="single" w:sz="4" w:space="0" w:color="999999"/><w:insideH w:val="single" w:sz="4" w:space="0" w:color="999999"/><w:insideV w:val="single" w:sz="4" w:space="0" w:color="999999"/></w:tblBorders></w:tblPr>'
  $xml += '<w:tr>'
  foreach ($h in $Headers) {
    $xml += '<w:tc><w:tcPr><w:shd w:fill="EDEDED"/></w:tcPr><w:p><w:r><w:rPr><w:b/><w:sz w:val="18"/></w:rPr><w:t>' + (Escape-XmlText $h) + '</w:t></w:r></w:p></w:tc>'
  }
  $xml += '</w:tr>'

  foreach ($row in $Rows) {
    $xml += '<w:tr>'
    foreach ($h in $Headers) {
      $v = $row.$h
      $xml += '<w:tc><w:p><w:r><w:rPr><w:sz w:val="18"/></w:rPr><w:t xml:space="preserve">' + (Escape-XmlText ([string]$v)) + '</w:t></w:r></w:p></w:tc>'
    }
    $xml += '</w:tr>'
  }
  $xml += '</w:tbl>'
  return $xml
}

function Get-ImageSizeInInches {
  param(
    [string]$Path,
    [double]$MaxWidthInches = 6.5
  )
  Add-Type -AssemblyName System.Drawing
  $img = [System.Drawing.Image]::FromFile($Path)
  try {
    $width = $img.Width
    $height = $img.Height
  } finally {
    $img.Dispose()
  }
  $ratio = $height / $width
  return @{ Width = $MaxWidthInches; Height = $MaxWidthInches * $ratio }
}

$buildRoot = Join-Path $env:TEMP ("conbio_docx_build_" + [guid]::NewGuid().ToString("N"))
$null = New-Item -ItemType Directory -Path $buildRoot
$null = New-Item -ItemType Directory -Path (Join-Path $buildRoot "_rels")
$null = New-Item -ItemType Directory -Path (Join-Path $buildRoot "word")
$null = New-Item -ItemType Directory -Path (Join-Path $buildRoot "word\_rels")
$null = New-Item -ItemType Directory -Path (Join-Path $buildRoot "word\media")
$null = New-Item -ItemType Directory -Path (Join-Path $buildRoot "docProps")

$body = New-Object System.Collections.Generic.List[string]
$body.Add((New-Paragraph "Class-specific spatiotemporal responses of mountain gazelles to predator context and land use" "Title"))
$body.Add((New-Paragraph "Raw Methods and Results package for Conservation Biology-style manuscript drafting" "Caption"))

$lines = Get-Content -LiteralPath $sourcePath
$skipSuggested = $false
$paragraphBuffer = New-Object System.Collections.Generic.List[string]

function Flush-ParagraphBuffer {
  if ($paragraphBuffer.Count -gt 0) {
    $text = ($paragraphBuffer -join " ").Trim()
    if ($text.Length -gt 0) {
      $body.Add((New-Paragraph $text "Normal"))
    }
    $paragraphBuffer.Clear()
  }
}

foreach ($line in $lines) {
  if ($line -like "## Suggested Figure and Table Package*") {
    Flush-ParagraphBuffer
    $skipSuggested = $true
    continue
  }
  if ($skipSuggested) { continue }

  $trim = $line.Trim()
  if ($trim.Length -eq 0) {
    Flush-ParagraphBuffer
    continue
  }
  if ($trim -match "^#\s+(.+)$") {
    Flush-ParagraphBuffer
    $body.Add((New-Paragraph $matches[1] "Heading1"))
  } elseif ($trim -match "^##\s+(.+)$") {
    Flush-ParagraphBuffer
    $body.Add((New-Paragraph $matches[1] "Heading1"))
  } elseif ($trim -match "^###\s+(.+)$") {
    Flush-ParagraphBuffer
    $body.Add((New-Paragraph $matches[1] "Heading2"))
  } elseif ($trim -match "^-\s+(.+)$") {
    Flush-ParagraphBuffer
    $body.Add((New-Paragraph $matches[1] "Bullet"))
  } else {
    $paragraphBuffer.Add($trim)
  }
}
Flush-ParagraphBuffer

$body.Add((New-PageBreak))
$body.Add((New-Paragraph "Figures" "Heading1"))

$figures = @(
  @{
    File = "publication_spabundance_wolf_jackal_relative_interaction_curves.png"
    Caption = "Figure 1. Class-specific relative abundance across wolf-jackal activity contexts. Curves show posterior medians and 95% credible intervals from separate negative-binomial spAbundance N-mixture models. Relative abundance is normalized to mean predator activity because abundance and detection intercepts did not converge."
  },
  @{
    File = "publication_fawn_female_beta_binomial_coefficients.png"
    Caption = "Figure 2. Predictors of fawn representation relative to adult females. Points show odds ratios and horizontal lines show 95% Wald confidence intervals from the primary season-adjusted additive beta-binomial mixed model. Continuous predictors were standardized."
  },
  @{
    File = "publication_predator_context_diel_by_class_median_combined.png"
    Caption = "Figure 3. Class-specific gazelle diel activity across four camera contexts defined by median wolf and jackal activity indices. Curves are normalized within class and context and therefore represent timing rather than event rate or abundance. The high-jackal/low-wolf fawn group was data deficient."
  },
  @{
    File = "publication_overall_gazelle_diel_activity.png"
    Caption = "Figure S1. Overall apparent-solar-time activity distributions of fawns, adult females, and adult males based on 30-min independent events. Curves are normalized activity densities."
  },
  @{
    File = "publication_spabundance_wolf_jackal_relative_heatmap.png"
    Caption = "Figure S2. Relative abundance across a 3 x 3 wolf-jackal interaction grid. Values are posterior median expected abundance divided by expected abundance at mean wolf and jackal activity."
  },
  @{
    File = "publication_predator_context_diel_pairwise_overlap_median.png"
    Caption = "Figure S3. Within-class temporal overlap among wolf-jackal contexts. Points show pairwise overlap coefficients and horizontal lines show 95% camera-cluster bootstrap intervals. Lower overlap indicates a larger difference in diel timing."
  }
)

$relationships = New-Object System.Collections.Generic.List[string]
$relationships.Add('<Relationship Id="rIdStyles" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>')
$relCounter = 1

foreach ($fig in $figures) {
  $fileName = $fig.File
  $imagePath = Join-Path $workspace $fileName
  if (-not (Test-Path -LiteralPath $imagePath)) {
    $body.Add((New-Paragraph ("Missing figure file: " + $fileName) "Caption"))
    continue
  }
  $relCounter += 1
  $relId = "rId" + $relCounter
  $targetName = "image" + $relCounter + ".png"
  Copy-Item -LiteralPath $imagePath -Destination (Join-Path $buildRoot ("word\media\" + $targetName))
  $relationships.Add('<Relationship Id="' + $relId + '" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/' + $targetName + '"/>')
  $size = Get-ImageSizeInInches -Path $imagePath -MaxWidthInches 6.5
  $body.Add((New-ImageParagraph -RelId $relId -Name $fileName -WidthInches $size.Width -HeightInches $size.Height))
  $body.Add((New-Paragraph $fig.Caption "Caption"))
}

$body.Add((New-PageBreak))
$body.Add((New-Paragraph "Key Tables" "Heading1"))

$spRowsRaw = Import-Csv -LiteralPath (Join-Path $workspace "publication_primary_spabundance_coefficients.csv")
$spRows = @()
foreach ($r in $spRowsRaw) {
  $spRows += [PSCustomObject]@{
    Class = $r.species
    Predictor = $r.predictor
    "Posterior mean [95% CrI]" = ("{0:0.00} [{1:0.00}, {2:0.00}]" -f [double]$r.estimate, [double]$r.lwr, [double]$r.upr)
    Rhat = ("{0:0.000}" -f [double]$r.rhat)
    ESS = ("{0:0}" -f [double]$r.ess)
    Supported = if ($r.credible_interval_excludes_zero -eq "TRUE") { "Yes" } else { "No" }
  }
}
$body.Add((New-Paragraph "Table 1. Class-specific wolf-jackal coefficients from the primary negative-binomial spAbundance N-mixture models." "Caption"))
$body.Add((New-WordTable -Headers @("Class", "Predictor", "Posterior mean [95% CrI]", "Rhat", "ESS", "Supported") -Rows $spRows))

$corRowsRaw = Import-Csv -LiteralPath (Join-Path $workspace "publication_focused_lfmsnmix_class_correlations.csv")
$corRows = @()
foreach ($r in $corRowsRaw) {
  $corRows += [PSCustomObject]@{
    "Class pair" = $r.class_pair
    "Residual correlation" = ("{0:0.00}" -f [double]$r.posterior_mean)
    "95% CrI" = ("{0:0.00}-{1:0.00}" -f [double]$r.lower_95, [double]$r.upper_95)
  }
}
$body.Add((New-Paragraph "Table 2. Residual correlations among demographic-class latent relative-abundance processes from the focused two-factor lfMsNMix model." "Caption"))
$body.Add((New-WordTable -Headers @("Class pair", "Residual correlation", "95% CrI") -Rows $corRows))

$betaRowsRaw = Import-Csv -LiteralPath (Join-Path $workspace "publication_fawn_female_beta_binomial_coefficients.csv")
$betaRows = @()
foreach ($r in $betaRowsRaw) {
  $pText = if ([double]$r.p_value -lt 0.001) { "<0.001" } else { "{0:0.###}" -f [double]$r.p_value }
  $betaRows += [PSCustomObject]@{
    Predictor = $r.predictor
    "Odds ratio [95% CI]" = ("{0:0.00} [{1:0.00}, {2:0.00}]" -f [double]$r.odds_ratio, [double]$r.conf_low, [double]$r.conf_high)
    P = $pText
  }
}
$body.Add((New-Paragraph "Table 3. Focal predictors of fawn representation relative to adult females from the primary beta-binomial mixed model." "Caption"))
$body.Add((New-WordTable -Headers @("Predictor", "Odds ratio [95% CI]", "P") -Rows $betaRows))

$dielRowsRaw = Import-Csv -LiteralPath (Join-Path $workspace "publication_authoritative_diel_context_comparisons.csv")
$dielRows = @()
foreach ($r in $dielRowsRaw) {
  $pText = if ($r.formal_test_eligible -eq "TRUE") {
    "{0:0.###}" -f [double]$r.compareCkern_p_fdr_within_class
  } else {
    "Data deficient"
  }
  $ciText = if ([string]::IsNullOrWhiteSpace($r.conf_low)) {
    "Not estimated"
  } else {
    ("{0:0.00}-{1:0.00}" -f [double]$r.conf_low, [double]$r.conf_high)
  }
  $dielRows += [PSCustomObject]@{
    Class = $r.class
    Contrast = $r.context_1 + " vs " + $r.context_2
    Estimator = $r.estimator
    Overlap = ("{0:0.00}" -f [double]$r.overlap)
    "95% bootstrap CI" = $ciText
    "FDR P" = $pText
  }
}
$body.Add((New-Paragraph "Table 4. Within-class diel overlap among median-split wolf-jackal camera contexts." "Caption"))
$body.Add((New-WordTable -Headers @("Class", "Contrast", "Estimator", "Overlap", "95% bootstrap CI", "FDR P") -Rows $dielRows))

$unmarkedRowsRaw = Import-Csv -LiteralPath (Join-Path $workspace "publication_unmarked_sensitivity_coefficients.csv")
$unmarkedRows = @()
foreach ($r in $unmarkedRowsRaw) {
  $pText = if ([double]$r.p_value -lt 0.001) { "<0.001" } else { "{0:0.###}" -f [double]$r.p_value }
  $unmarkedRows += [PSCustomObject]@{
    Class = $r.species
    Predictor = $r.predictor
    "Estimate [95% CI]" = ("{0:0.00} [{1:0.00}, {2:0.00}]" -f [double]$r.estimate, [double]$r.conf_low, [double]$r.conf_high)
    P = $pText
  }
}
$body.Add((New-Paragraph "Table S1. Matched unmarked negative-binomial pcount sensitivity estimates." "Caption"))
$body.Add((New-WordTable -Headers @("Class", "Predictor", "Estimate [95% CI]", "P") -Rows $unmarkedRows))

$documentXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:wpc="http://schemas.microsoft.com/office/word/2010/wordprocessingCanvas"
  xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006"
  xmlns:o="urn:schemas-microsoft-com:office:office"
  xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
  xmlns:m="http://schemas.openxmlformats.org/officeDocument/2006/math"
  xmlns:v="urn:schemas-microsoft-com:vml"
  xmlns:wp14="http://schemas.microsoft.com/office/word/2010/wordprocessingDrawing"
  xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
  xmlns:w10="urn:schemas-microsoft-com:office:word"
  xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"
  xmlns:w14="http://schemas.microsoft.com/office/word/2010/wordml"
  xmlns:wpg="http://schemas.microsoft.com/office/word/2010/wordprocessingGroup"
  xmlns:wpi="http://schemas.microsoft.com/office/word/2010/wordprocessingInk"
  xmlns:wne="http://schemas.microsoft.com/office/word/2006/wordml"
  xmlns:wps="http://schemas.microsoft.com/office/word/2010/wordprocessingShape"
  mc:Ignorable="w14 wp14">
  <w:body>
    $($body -join "`n")
    <w:sectPr>
      <w:pgSz w:w="12240" w:h="15840"/>
      <w:pgMar w:top="1440" w:right="1080" w:bottom="1440" w:left="1080" w:header="720" w:footer="720" w:gutter="0"/>
    </w:sectPr>
  </w:body>
</w:document>
"@

$stylesXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:style w:type="paragraph" w:default="1" w:styleId="Normal">
    <w:name w:val="Normal"/>
    <w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="22"/></w:rPr>
  </w:style>
</w:styles>
"@

$relsXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>
"@

$documentRelsXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  $($relationships -join "`n  ")
</Relationships>
"@

$contentTypesXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Default Extension="png" ContentType="image/png"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
  <Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>
"@

$created = (Get-Date).ToUniversalTime().ToString("s") + "Z"
$coreXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties"
  xmlns:dc="http://purl.org/dc/elements/1.1/"
  xmlns:dcterms="http://purl.org/dc/terms/"
  xmlns:dcmitype="http://purl.org/dc/dcmitype/"
  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:title>Golan Gazelles Conservation Biology Methods and Results</dc:title>
  <dc:creator>Codex analysis workflow</dc:creator>
  <cp:lastModifiedBy>Codex analysis workflow</cp:lastModifiedBy>
  <dcterms:created xsi:type="dcterms:W3CDTF">$created</dcterms:created>
  <dcterms:modified xsi:type="dcterms:W3CDTF">$created</dcterms:modified>
</cp:coreProperties>
"@

$appXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"
  xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>Codex OOXML builder</Application>
</Properties>
"@

Set-Content -LiteralPath (Join-Path $buildRoot "[Content_Types].xml") -Value $contentTypesXml -Encoding UTF8
Set-Content -LiteralPath (Join-Path $buildRoot "_rels\.rels") -Value $relsXml -Encoding UTF8
Set-Content -LiteralPath (Join-Path $buildRoot "word\document.xml") -Value $documentXml -Encoding UTF8
Set-Content -LiteralPath (Join-Path $buildRoot "word\styles.xml") -Value $stylesXml -Encoding UTF8
Set-Content -LiteralPath (Join-Path $buildRoot "word\_rels\document.xml.rels") -Value $documentRelsXml -Encoding UTF8
Set-Content -LiteralPath (Join-Path $buildRoot "docProps\core.xml") -Value $coreXml -Encoding UTF8
Set-Content -LiteralPath (Join-Path $buildRoot "docProps\app.xml") -Value $appXml -Encoding UTF8

if (Test-Path -LiteralPath $outputPath) {
  Remove-Item -LiteralPath $outputPath -Force
}

$fileStream = [System.IO.File]::Open($outputPath, [System.IO.FileMode]::Create)
try {
  $archive = New-Object System.IO.Compression.ZipArchive($fileStream, [System.IO.Compression.ZipArchiveMode]::Create)
  try {
    foreach ($file in Get-ChildItem -LiteralPath $buildRoot -Recurse -File) {
      $entryName = $file.FullName.Substring($buildRoot.Length + 1).Replace('\', '/')
      $entry = $archive.CreateEntry($entryName, [System.IO.Compression.CompressionLevel]::Optimal)
      $entryStream = $entry.Open()
      $sourceStream = [System.IO.File]::OpenRead($file.FullName)
      try {
        $sourceStream.CopyTo($entryStream)
      } finally {
        $sourceStream.Dispose()
        $entryStream.Dispose()
      }
    }
  } finally {
    $archive.Dispose()
  }
} finally {
  $fileStream.Dispose()
}

Write-Host "Wrote $outputPath"
