# API Reference

Documentation for `Endas.jl`'s public interface.

## Contents

```@contents
Pages = ["reference.md"]
Depth = 3
```

## Modules

```@autodocs
Modules = [Endas]
Order   = [:module]
```

## Types

```@index
Modules = [Endas]
Order   = [:type]
```

### Abstract types

```@autodocs
Modules = [Endas]
Order   = [:type]
Filter  = t -> isabstracttype(t)
```

### Observation operators

```@autodocs
Modules = [Endas]
Order   = [:type]
Filter  = t -> t <: ObservationOperator && !isabstracttype(t)
```

### Covariance operators

```@autodocs
Modules = [Endas]
Order   = [:type]
Filter  = t -> t <: CovarianceOperator && !isabstracttype(t)
```

### Filters and smoothers

```@autodocs
Modules = [Endas]
Order   = [:type]
Filter  = t -> t <: AbstractFilter && !isabstracttype(t)
```

### Models

```@autodocs
Modules = [Endas]
Order   = [:type]
Filter  = t -> t <: AbstractModel && !isabstracttype(t)
```


## Functions

```@index
Modules = [Endas]
Order   = [:function]
```

```@autodocs
Modules = [Endas]
Order   = [:function]
```
