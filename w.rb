def find_one_using_hash_map(array)
  map = {}
  dup = nil
  array.each do |v|
    map[v] = (map[v] || 0 ) + 1

    if map[v] > 1
      dup = v
      break
    end
  end

  return dup
end

f1 = File.read('6_1.txt').to_s.split("\n").map do |s|
	s.partition('=').last
end
f2 = File.read('6_4.txt').to_s.split("\n").map do |s|
	s.partition('=').last
end

puts f1 - f2
# puts f1.uniq.count
# puts f2.uniq.count
# puts f2.uniq.count
puts find_one_using_hash_map(f1)
puts find_one_using_hash_map(f2)
puts "-----"
puts f2 - f1
# puts f1
