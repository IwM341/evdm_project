#pragma once

#include <variant>
#include <vector>
#include <memory>
#include <cmath>
#include <array>

namespace evdm {

	struct mutex_empty_t {
		inline void lock() {}
		inline void unlock() {}
	};
	template<typename Fp_t, typename T, typename U = bool, size_t max_lvl = 30>
	struct Quadtree {
		typedef uint_least64_t u64_t;
		struct Node {
		private:
			//(0,0)=(left down) (0,1)=(left up), (1,0) = (right down), (1,1) = (right up)
			std::array<std::unique_ptr<Node>, 4> children = { nullptr,nullptr,nullptr,nullptr }; // NW, NE, SW, SE
			Node* _parent = nullptr;
			u64_t _level = 0;
			u64_t cx0, cx1, cy0, cy1;
			Fp_t _x0, _x1, _y0, _y1;
			std::array<std::shared_ptr<T>, 4> corner_values; // NW, NE, SW, SE corners
			bool _is_leaf;

			friend struct Quadtree;
			Node(
				u64_t cx0, u64_t cx1, u64_t cy0, u64_t cy1,
				Quadtree& m_tree, Node* m_parent = nullptr, u64_t m_level = 0)
				: cx0(cx0), cx1(cx1), cy0(cy0), cy1(cy1),
				_x0(m_tree.x_val(cx0)), _x1(m_tree.x_val(cx1)),
				_y0(m_tree.y_val(cy0)), _y1(m_tree.y_val(cy1)),
				_parent(m_parent), _level(m_level), _is_leaf(true) {}


		public:
#define DECLARE_CONST_PROPERTY(Type,prop) \
		inline Type prop() const {return _##prop;}
			DECLARE_CONST_PROPERTY(Fp_t, x0)
				DECLARE_CONST_PROPERTY(Fp_t, x1)
				DECLARE_CONST_PROPERTY(Fp_t, y0)
				DECLARE_CONST_PROPERTY(Fp_t, y1)
#undef DECLARE_CONST_PROPERTY
				U value;
			size_t level()const {
				return _level;
			}
			bool is_leaf()const {
				return _is_leaf;
			}
			Node& child(size_t i) {
				if (_is_leaf) {
					throw std::runtime_error("leaf has no child");
				}
				return *children[i];
			}
			Node const& child(size_t i)const {
				if (_is_leaf) {
					throw std::runtime_error("leaf has no child");
				}
				return *children[i];
			}
			Node& parent() {
				if (_parent == nullptr) {
					throw std::runtime_error("root has no parent");
				}
				return *_parent;
			}
			bool contain(Fp_t x, Fp_t y)const {
				return _x0 <= x && x <= _x1 && _y0 <= y && y <= _y1;
			}
			std::pair<Fp_t, Fp_t> coeffs(Fp_t x, Fp_t y)const {
				return { (x1() - x) / (x1() - x0()),(y1() - y) / (y1() - y0()) };
			}
			const T& operator ()(size_t i) const {
				return *corner_values[i];
			}
			T& operator ()(size_t i) {
				return *corner_values[i];
			}
			const T& operator ()(size_t i, size_t j) const {
				return *corner_values[2 * i + j];
			}
			T& operator ()(size_t i, size_t j) {
				return *corner_values[2 * i + j];
			}

		};

		std::unique_ptr<Node> m_root;
		std::map<u64_t, std::shared_ptr<T>> point_registry;

		Fp_t X0, X1, Y0, Y1;
		u64_t max_level_dyn = 0;
		constexpr static u64_t _1L = ((u64_t)1);
		constexpr static u64_t _3L = ((u64_t)3);
		constexpr static u64_t max_c = _1L << max_lvl;
		constexpr static u64_t max_c_inner = (_1L << max_lvl) - 1;

		Fp_t dX_inv() const {
			return 1 / (X1 - X0);
		}
		Fp_t dY_inv() const {
			return 1 / (Y1 - Y0);
		}
		// Хеш функция для точки
		u64_t point_coord(Fp_t x, Fp_t y) const {
			u64_t hx = x_coord(x);
			u64_t hy = y_coord(y);
			return (hx << (max_lvl + 1)) | hy;
		}

		u64_t x_coord(Fp_t x)const {
			auto max_c_d = max_c;
			return static_cast<u64_t>(((x - X0) * dX_inv()) * max_c);
		}
		u64_t y_coord(Fp_t y)const {
			auto max_c_d = max_c;
			return static_cast<uint_least64_t>(((y - Y0) * dY_inv()) * max_c);
		}

		Fp_t x_val(uint64_t xc)const {
			Fp_t a = (Fp_t(xc)) / max_c;
			return X0 + a * (X1 - X0);
		}
		Fp_t y_val(u64_t yc)const {
			Fp_t b = (Fp_t(yc)) / max_c;
			return Y0 + b * (Y1 - Y0);
		}
		u64_t get_key(u64_t xc, u64_t yc) {
			return (xc << (max_lvl + 1)) | yc;
		}
		std::pair<u64_t, u64_t> xcyc(u64_t key) {
			return { key >> (max_lvl + 1),key & max_c_inner };
		}

		// Получение или создание значения для точки
		template <typename Value_Gen_t, typename Mutex_t = mutex_empty_t >
		std::shared_ptr<T> get_or_create_value(
			u64_t xc, u64_t yc, 
			Value_Gen_t&& v_creator,
			Mutex_t &&m_mutex = mutex_empty_t{}) {
			auto m_point = get_key(xc, yc);
			{
				std::lock_guard<Mutex_t> m_guard(m_mutex);
				if (auto it = point_registry.find(m_point); it != point_registry.end()) {
					return it->second;
				}
			}
			auto value = std::make_shared<T>(v_creator());
			{
				std::lock_guard<Mutex_t> m_guard(m_mutex);
				point_registry[m_point] = value;
			}

			return value;
		}



		// Наследование значений углов от родителя
		template <typename Fun_t,
			typename Mutex_t = mutex_empty_t >
		void inherit_corner_values(Node& node, Fun_t&& F, 
			Mutex_t&& m_mutex = mutex_empty_t{}) {

			node.corner_values[0] = get_or_create_value(node.cx0, node.cy0,
				[x = node.x0(), y = node.y0(), &F]() {return F(x, y); }, m_mutex);

			node.corner_values[1] = get_or_create_value(node.cx0, node.cy1,
				[x = node.x0(), y = node.y1(), &F]() {return F(x, y); }, m_mutex);

			node.corner_values[2] = get_or_create_value(node.cx1, node.cy0,
				[x = node.x1(), y = node.y0(), &F]() {return F(x, y); }, m_mutex);

			node.corner_values[3] = get_or_create_value(node.cx1, node.cy1,
				[x = node.x1(), y = node.y1(), &F]() {return F(x, y); }, m_mutex);

		}


		std::unique_ptr<Node> create_node(
			u64_t x0, u64_t x1, u64_t y0, u64_t y1,
			Node* parent = nullptr, int level = 0)
		{
			return std::make_unique<Node>(Node(x0, x1, y0, y1, *this, parent, level));
		}

		const Quadtree& cthis()const {
			return *this;
		}


	public:
		template <typename FillFunc>
		Quadtree(FillFunc&& F, Fp_t X0 = 0, Fp_t X1 = 1, Fp_t Y0 = 0, Fp_t Y1 = 1)
			:X0(X0), X1(X1), Y0(Y0), Y1(Y1)
		{
			static_assert(max_lvl <= 31, "max_lvl should be less 31");
			if (X1 <= X0 || Y1 <= Y0) {
				throw std::runtime_error("X1 <= X0 || Y1 <= Y0");
			}
			m_root = create_node(0, max_c, 0, max_c);
			inherit_corner_values(*m_root, F);
		}
		Node& root() {
			return *m_root;
		}
		Node const& root()const {
			return *m_root;
		}
		template <typename Value_Gen_t>
		void refine(Node& node, size_t num_levels, Value_Gen_t&& F) {
			if (num_levels == 0) {
				return;
			}
			else {
				if (node.is_leaf()) {
					split_node(node, F);
					for (auto& ch : node.children) {
						refine(*ch, num_levels - 1, F);
					}
				}
				else {
					for (auto& ch : node.children) {
						refine(*ch, num_levels, F);
					}
				}
			}
		}

		template <typename Value_Gen_t>
		void refine_to_level(Node& node, size_t max_level, Value_Gen_t&& F) {
			if (max_level == 0) {
				return;
			}
			else {
				if (node.is_leaf()) {
					split_node(node, F);
					for (auto& ch : node.children) {
						refine(*ch, max_level - 1, F);
					}
				}
				else {
					for (auto& ch : node.children) {
						refine(*ch, max_level - 1, F);
					}
				}
			}
		}

		// Разделение узла
		template <typename Value_Gen_t>
		void split_node(Node& node, Value_Gen_t&& value_func) {
			if (!node.is_leaf()) return;
			if (node.level() >= max_lvl) {
				throw std::runtime_error("max level reached");
			}



			Fp_t cx_n = (node.cx0 + node.cx1) / 2;
			Fp_t cy_n = (node.cy0 + node.cy1) / 2;

			// Создаем дочерние узлы
			u64_t new_level = node.level() + 1;
			auto n00 = create_node(node.cx0, cx_n, node.cy0, cy_n, &node, new_level);
			auto n01 = create_node(node.cx0, cx_n, cy_n, node.cy1, &node, new_level);
			auto n10 = create_node(cx_n, node.cx1, node.cy0, cy_n, &node, new_level);
			auto n11 = create_node(cx_n, node.cx1, cy_n, node.cy1, &node, new_level);
			node.children[0] = std::move(n00);
			node.children[1] = std::move(n01);
			node.children[2] = std::move(n10);
			node.children[3] = std::move(n11);
			node._is_leaf = false;


			// Наследуем значения углов от родителя
			inherit_corner_values(*node.children[0], value_func);
			inherit_corner_values(*node.children[1], value_func);
			inherit_corner_values(*node.children[2], value_func);
			inherit_corner_values(*node.children[3], value_func);

			max_level_dyn = std::max(max_level_dyn, new_level);

		}


		// Быстрый поиск узла содержащего точку
		inline Node const& find_node_coord(u64_t cx, u64_t cy) const {
			Node* tmp_node = m_root.get();
			while (!tmp_node->is_leaf()) {
				u64_t ix = 1 & (cx >> (max_lvl - tmp_node->_level - 1));
				u64_t iy = 1 & (cy >> (max_lvl - tmp_node->_level - 1));
				tmp_node = &tmp_node->child(ix << 1 | iy);
			}
			return *tmp_node;
		}
		inline Node& find_node_coord(u64_t cx, u64_t cy) {
			return const_cast<Node&>(cthis().find_node_coord(cx, cy));
		}

		Node const& find_node(Fp_t x, Fp_t y) const {
			x = std::clamp(x, X0, X1);
			y = std::clamp(y, Y0, Y1);

			u64_t cx = std::min(x_coord(x), max_c_inner);
			u64_t cy = std::min(y_coord(y), max_c_inner);
			return find_node_coord(cx, cy);
		}
		inline Node& find_node(Fp_t x, Fp_t y) {
			return const_cast<Node&>(cthis().find_node(x, y));
		}

		bool contain(Fp_t x, Fp_t y) const {
			return root().contain(x, y);
		}

		template <typename Node_t>
		struct leaf_iterator_t {
			Node_t* m_node;

			leaf_iterator_t(Node_t* m_node) :m_node(m_node) {}
			inline operator bool()const {
				return m_node != nullptr;
			}
			inline bool operator ==(leaf_iterator_t const& other)const {
				return m_node == other.m_node;
			}
			inline leaf_iterator_t& operator++() {
				u64_t morton_code = detail::morton_encode(
					m_node->cx0, m_node->cy0
				);
				/*
				println("tmpcode");
				println(std::bitset<8>(morton_code));
				println(std::bitset<4>(m_node->cx0));
				println(std::bitset<4>(m_node->cy0));
				println("tmp codes");
				std::cout << std::bitset<31>(m_node->cx0) <<
					"\n" << std::bitset<31>(m_node->cy0) << std::endl;
				std::cout << std::bitset<62>(morton_code)<< std::endl;

				u64_t morton_code_new = detail::morton_encode(
							m_node->_parent->children[1]->cx0,
							m_node->_parent->children[1]->cy0
						);
				println("desire code");
				std::cout << std::bitset<31>(
					m_node->_parent->children[1]->cx0) << std::endl;

				std::cout << std::bitset<31>(
					m_node->_parent->children[1]->cy0) << std::endl;

				std::cout << std::bitset<62>(morton_code_new)<< std::endl;
				*/
				u64_t new_code = morton_code + (_1L << (2 * (max_lvl - m_node->_level)));

				//println("new code");
				//println(std::bitset<8>(new_code));
				/*println("new code");
				std::cout << std::bitset<62>(new_code)<< std::endl;
				*/
				//assert(cy == );
				Node_t* common_parent = find_common_parent(
					m_node, morton_code, new_code,
					(_3L) << (2 * (max_lvl - m_node->_level))
				);
				/*println("mask");
				std::cout << std::bitset<62>(
					(_3L) << (2*(max_lvl - m_node->_level))
				)<< std::endl;*/
				if (common_parent == nullptr) {
					m_node = nullptr;
					return *this;
				}

				m_node = find_leaf(common_parent, new_code,
					2 * (max_lvl - common_parent->_level - 1));

				//assert(m_node == &find_node_coord(cx,cy));
				return *this;
			};
			inline Node_t& operator *() {
				return m_node;
			}
			inline Node_t* operator ->() {
				return m_node;
			}

		private:
			friend struct Quadtree;
			static inline Node_t* find_common_parent(
				Node_t* tmp_node, u64_t zcode, u64_t tmp_code, u64_t level_mask) {
				while ((zcode & level_mask) != (tmp_code & level_mask)) {
					tmp_node = tmp_node->_parent;
					level_mask = level_mask << 2;
				}
				return tmp_node;
			}
			static inline Node_t* find_leaf(
				Node_t* parent, u64_t m_code, u64_t m_level) {
				while (!parent->is_leaf()) {
					size_t index = (m_code >> m_level) & _3L;
					parent = &(parent->child(index));
					m_level -= 2;
				}
				return parent;
			}
		};

		leaf_iterator_t<Node> begin() {
			return { &find_node_coord(0,0) };
		}

		leaf_iterator_t<Node> end() {
			return { nullptr };
		}

	};


};