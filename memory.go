package sparse

import "github.com/Konstantin8105/mms"

var (
	floats = mms.Float64sCache{}
	ints   = mms.IntsCache{}
)

func init() {
	mms.Debug = false
}
